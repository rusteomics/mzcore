//! LZW decompression for classic Unix .Z files
//!
//! Based on the classic Unix "uncompress" / "zcat" implementation (ncompress).
//! Supports 9-16 bit codes and block mode (CLEAR code).
//!
//! The .Z format packs codes in blocks of `n_bits` bytes (8 codes per block).
//! When a CLEAR code is emitted or the code width increases, the compressor
//! pads the remainder of the current block. The decompressor must therefore
//! discard the padding by reading codes in aligned blocks, matching the
//! original ncompress `getcode()` behavior.

use crate::lzw::ArchiveError;

/// Maximum bits for codes (9-16)
const BITS: usize = 16;
/// Dictionary table size
const HSIZE: usize = 69001;
/// Initial code width
const INIT_BITS: usize = 9;

/// Max code for a given bit width
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
        let maxmaxcode: usize = 1 << max_bits;

        // Block-based code reading state (matches ncompress getcode() behavior).
        // Codes are packed in blocks of n_bits bytes (= 8 codes). When n_bits
        // changes or a CLEAR code is encountered, the remaining bits in the
        // current block are padding and must be discarded.
        let mut block_buf = [0u8; BITS];
        let mut block_bit_offset: usize = 0;
        let mut block_bit_size: usize = 0;
        let mut input_pos: usize = 0;
        let mut clear_flag = false;

        // LZW state
        let mut n_bits: usize = INIT_BITS;
        let mut maxcode_val = maxcode(n_bits);
        let mut free_ent: usize = if block_mode { 257 } else { 256 };

        // Dictionary tables
        let mut tab_prefix: Vec<u16> = vec![0; HSIZE];
        let mut tab_suffix: Vec<u8> = vec![0; HSIZE];

        // Output buffer
        let mut output: Vec<u8> = Vec::new();

        // Stack for decoding (strings are built in reverse)
        let mut stack: Vec<u8> = Vec::with_capacity(8192);

        // -- getcode! macro: reads the next code from the input using block-
        //    aligned reading, exactly as ncompress does. When clear_flag is set,
        //    free_ent > maxcode_val, or the current block is exhausted, a fresh
        //    block of n_bits bytes is read from the input. --
        macro_rules! getcode {
            () => {{
                // Check whether we need to read a new block
                if clear_flag || block_bit_offset >= block_bit_size || free_ent > maxcode_val {
                    // Increase code width if the dictionary outgrew it
                    if free_ent > maxcode_val {
                        n_bits += 1;
                        if n_bits == max_bits {
                            maxcode_val = maxmaxcode;
                        } else {
                            maxcode_val = maxcode(n_bits);
                        }
                    }
                    // After CLEAR, reset code width
                    if clear_flag {
                        n_bits = INIT_BITS;
                        maxcode_val = maxcode(n_bits);
                        clear_flag = false;
                    }
                    // Read a fresh block of n_bits bytes
                    let avail = input.len().saturating_sub(input_pos);
                    let to_read = n_bits.min(avail);
                    if to_read == 0 {
                        None::<usize>
                    } else {
                        block_buf[..to_read]
                            .copy_from_slice(&input[input_pos..input_pos + to_read]);
                        input_pos += to_read;
                        block_bit_offset = 0;
                        let bit_count = (to_read << 3) as isize - (n_bits as isize - 1);
                        if bit_count <= 0 {
                            None::<usize>
                        } else {
                            block_bit_size = bit_count as usize;
                            // fall through to extraction below
                            let bp = block_bit_offset >> 3;
                            let r_off = block_bit_offset & 7;
                            let mut code = (block_buf[bp] as usize) >> r_off;
                            let mut bits_got = 8 - r_off;
                            if bits_got < n_bits {
                                code |= (block_buf[bp + 1] as usize) << bits_got;
                                bits_got += 8;
                            }
                            if bits_got < n_bits {
                                code |= (block_buf[bp + 2] as usize) << bits_got;
                            }
                            code &= (1usize << n_bits) - 1;
                            block_bit_offset += n_bits;
                            Some(code)
                        }
                    }
                } else {
                    // Extract code from current block
                    let bp = block_bit_offset >> 3;
                    let r_off = block_bit_offset & 7;
                    let mut code = (block_buf[bp] as usize) >> r_off;
                    let mut bits_got = 8 - r_off;
                    if bits_got < n_bits {
                        code |= (block_buf[bp + 1] as usize) << bits_got;
                        bits_got += 8;
                    }
                    if bits_got < n_bits {
                        code |= (block_buf[bp + 2] as usize) << bits_got;
                    }
                    code &= (1usize << n_bits) - 1;
                    block_bit_offset += n_bits;
                    Some(code)
                }
            }};
        }

        // Read first code (must be a literal byte)
        let first = match getcode!() {
            Some(c) => c,
            None => return Ok(output),
        };

        if first > 255 {
            return Err(ArchiveError::DecompressionFailed {
                entry: String::new(),
                reason: "First code > 255".to_string(),
            });
        }

        let mut finchar = first as u8;
        let mut oldcode = first;
        output.push(finchar);

        // Main decompression loop
        loop {
            let incode = match getcode!() {
                Some(c) => c,
                None => break,
            };

            // Handle CLEAR code in block mode
            if incode == 256 && block_mode {
                // Reset dictionary
                tab_prefix.fill(0);
                free_ent = 257;
                clear_flag = true;
                // After CLEAR, the stream restarts with a literal code that
                // initializes oldcode/finchar. This call sees clear_flag,
                // realigns to the next block, and resets code width.
                let next = match getcode!() {
                    Some(c) => c,
                    None => break,
                };

                if next > 255 {
                    return Err(ArchiveError::DecompressionFailed {
                        entry: String::new(),
                        reason: format!("Code {} after CLEAR is not a literal", next),
                    });
                }

                finchar = next as u8;
                oldcode = next;
                output.push(finchar);
                continue;
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
                stack.push(finchar);
                code = oldcode;
            }

            // Decode the string by following the chain
            while code >= 256 {
                stack.push(tab_suffix[code]);
                code = tab_prefix[code] as usize;
            }

            // First character of the string
            finchar = code as u8;

            // Output the decoded string (in correct order)
            output.push(finchar);
            while let Some(c) = stack.pop() {
                output.push(c);
            }

            // Add new entry to dictionary
            if free_ent < maxmaxcode {
                tab_prefix[free_ent] = oldcode as u16;
                tab_suffix[free_ent] = finchar;
                free_ent += 1;
            }

            oldcode = incode;
        }

        Ok(output)
    }
}
