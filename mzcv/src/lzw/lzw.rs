//! LZW decompression for classic Unix .Z files
//!
//! Based on the classic Unix "uncompress" / "zcat" implementation. (from ncompress)
//! Supports 9-16 bit codes and block mode (CLEAR code).

use crate::lzw::ArchiveError;

/// Maximum bits for codes (9-16)
const BITS: usize = 16;
/// Hash table size
const HSIZE: usize = 69001;

/// Macro to calculate max code for n bits
#[inline]
fn maxcode(n: usize) -> usize {
    (1 << n) - 1
}

/// LZW decompressor for .Z files
pub(crate) struct Lzw {
    max_bits: u8,
    block_mode: bool,
}

impl Lzw {
    /// Create a new LZW decompressor
    ///
    /// # Arguments
    /// * `max_bits` - Maximum bits per code (9-16), from .Z file header
    /// * `block_mode` - If true, handle CLEAR codes (code 256) for table reset
    pub(crate) fn new(max_bits: u8, block_mode: bool) -> Self {
        Lzw {
            max_bits,
            block_mode,
        }
    }

    /// Decompress LZW-compressed data
    pub(crate) fn decomp(&mut self, input: &[u8]) -> Result<Vec<u8>, ArchiveError> {
        let max_bits = self.max_bits as usize;
        let block_mode = self.block_mode;

        // Bit buffer for reading variable-width codes
        let mut buf: u64 = 0;
        let mut buflen: usize = 0;
        let mut input_pos: usize = 0;

        // Read next code from input stream
        let read_code = |n_bits: usize,
                         buf: &mut u64,
                         buflen: &mut usize,
                         input_pos: &mut usize|
         -> Option<usize> {
            while *buflen < n_bits {
                if *input_pos >= input.len() {
                    return None;
                }
                *buf |= (input[*input_pos] as u64) << *buflen;
                *input_pos += 1;
                *buflen += 8;
            }
            let code = (*buf & ((1u64 << n_bits) - 1)) as usize;
            *buf >>= n_bits;
            *buflen -= n_bits;
            Some(code)
        };

        // Dictionary tables
        let mut htab: Vec<u16> = vec![0; HSIZE]; // prefix table
        let mut codetab: Vec<u8> = vec![0; HSIZE]; // suffix table

        // Output buffer
        let mut output: Vec<u8> = Vec::new();

        // Stack for decoding (strings are built in reverse)
        let mut stack: Vec<u8> = Vec::with_capacity(HSIZE);

        // Initialize
        let mut n_bits = 9usize;
        let mut maxcode_val = maxcode(n_bits);
        let mut free_ent = if block_mode { 257 } else { 256 };

        // Read first code
        let oldcode = match read_code(n_bits, &mut buf, &mut buflen, &mut input_pos) {
            Some(c) => c,
            None => return Ok(output),
        };

        if oldcode > 255 {
            return Err(ArchiveError::DecompressionFailed {
                entry: String::new(),
                reason: "First code > 255".to_string(),
            });
        }

        let mut finchar = oldcode as u8;
        let mut oldcode = oldcode;
        output.push(finchar);

        // Main decompression loop
        while let Some(incode) = read_code(n_bits, &mut buf, &mut buflen, &mut input_pos) {
            // Handle CLEAR code in block mode
            if incode == 256 && block_mode {
                // Reset dictionary
                htab.fill(0);
                n_bits = 9;
                maxcode_val = maxcode(n_bits);

                // Read next code after clear
                match read_code(n_bits, &mut buf, &mut buflen, &mut input_pos) {
                    Some(c) => {
                        if c > 255 {
                            return Err(ArchiveError::DecompressionFailed {
                                entry: String::new(),
                                reason: "Code after CLEAR > 255".to_string(),
                            });
                        }
                        finchar = c as u8;
                        oldcode = c;
                        output.push(finchar);
                        free_ent = 257;
                        continue;
                    }
                    None => break,
                }
            }

            let mut code = incode;

            // Special case: code == free_ent (KwKwK)
            if code >= free_ent {
                if code > free_ent {
                    return Err(ArchiveError::DecompressionFailed {
                        entry: String::new(),
                        reason: format!("Invalid code {} > free_ent {}", code, free_ent),
                    });
                }
                // Push finchar and use oldcode
                stack.push(finchar);
                code = oldcode;
            }

            // Decode the string by following the chain
            while code >= 256 {
                stack.push(codetab[code]);
                code = htab[code] as usize;
            }

            // First character of the string
            finchar = code as u8;

            // Output the decoded string (in correct order)
            output.push(finchar);
            while let Some(c) = stack.pop() {
                output.push(c);
            }

            // Add new entry to dictionary
            if free_ent < (1 << BITS) {
                htab[free_ent] = oldcode as u16;
                codetab[free_ent] = finchar;
                free_ent += 1;

                // Increase code size if needed
                if free_ent > maxcode_val && n_bits < max_bits {
                    n_bits += 1;
                    maxcode_val = maxcode(n_bits);
                }
            }

            oldcode = incode;
        }

        Ok(output)
    }
}
