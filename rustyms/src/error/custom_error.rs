use super::Context;
use itertools::Itertools;
use serde::*;
use std::error;
use std::fmt;

/// An error. Stored as a pointer to a structure on the heap to prevent large sizes which could be
/// detrimental to performance for the happy path.
#[derive(Serialize, Deserialize, PartialEq, Clone, Eq, Hash)]
pub struct CustomError {
    content: Box<InnerError>,
}

#[derive(Serialize, Deserialize, PartialEq, Clone, Eq, Hash)]
struct InnerError {
    /// The level of the error, defining how it should be handled
    warning: bool,
    /// A short description of the error, generally used as title line
    short_description: String,
    /// A longer description of the error, presented below the context to give more information and helpful feedback
    long_description: String,
    /// Possible suggestion(s) for the indicated text
    suggestions: Vec<String>,
    /// Version if applicable
    version: String,
    /// The context, in the most general sense this produces output which leads the user to the right place in the code or file
    context: Context,
    /// Underlying errors
    underlying_errors: Vec<CustomError>,
}

#[expect(clippy::needless_pass_by_value)] // The impl ToString should be passed like this, otherwise &str gives errors
impl CustomError {
    /// Create a new `CustomError`.
    ///
    /// ## Arguments
    /// * `short_desc` - A short description of the error, generally used as title line.
    /// * `long_desc` -  A longer description of the error, presented below the context to give more information and helpful feedback.
    /// * `context` - The context, in the most general sense this produces output which leads the user to the right place in the code or file.
    pub fn error(short_desc: impl ToString, long_desc: impl ToString, context: Context) -> Self {
        Self {
            content: Box::new(InnerError {
                warning: false,
                short_description: short_desc.to_string(),
                long_description: long_desc.to_string(),
                suggestions: Vec::new(),
                version: String::new(),
                context,
                underlying_errors: Vec::new(),
            }),
        }
    }
    /// Create a new `CustomError`.
    ///
    /// ## Arguments
    /// * `short_desc` - A short description of the error, generally used as title line.
    /// * `long_desc` -  A longer description of the error, presented below the context to give more information and helpful feedback.
    /// * `context` - The context, in the most general sense this produces output which leads the user to the right place in the code or file.
    pub fn warning(short_desc: impl ToString, long_desc: impl ToString, context: Context) -> Self {
        Self {
            content: Box::new(InnerError {
                warning: true,
                short_description: short_desc.to_string(),
                long_description: long_desc.to_string(),
                suggestions: Vec::new(),
                version: String::new(),
                context,
                underlying_errors: Vec::new(),
            }),
        }
    }

    /// The level of the error
    pub const fn level(&self) -> &str {
        if self.content.warning {
            "warning"
        } else {
            "error"
        }
    }

    /// The suggestions
    pub fn suggestions(&self) -> &[String] {
        &self.content.suggestions
    }

    /// Tests if this errors is a warning
    pub const fn is_warning(&self) -> bool {
        self.content.warning
    }

    /// Gives the short description or title for this error
    pub fn short_description(&self) -> &str {
        &self.content.short_description
    }

    /// Gives the long description for this error
    pub fn long_description(&self) -> &str {
        &self.content.long_description
    }

    /// Create a copy of the error with a new long description
    #[must_use]
    pub fn with_long_description(&self, long_desc: impl ToString) -> Self {
        Self {
            content: Box::new(InnerError {
                long_description: long_desc.to_string(),
                ..(*self.content).clone()
            }),
        }
    }

    /// Create a copy of the error with the given suggestions
    #[must_use]
    pub fn with_suggestions(&self, suggestions: impl IntoIterator<Item = impl ToString>) -> Self {
        Self {
            content: Box::new(InnerError {
                suggestions: suggestions.into_iter().map(|s| s.to_string()).collect(),
                ..(*self.content).clone()
            }),
        }
    }

    /// Set the version of the underlying format
    #[must_use]
    pub fn with_version(self, version: impl ToString) -> Self {
        Self {
            content: Box::new(InnerError {
                version: version.to_string(),
                ..(*self.content)
            }),
        }
    }

    /// Create a copy of the error with a new context
    #[must_use]
    pub fn with_context(&self, context: Context) -> Self {
        Self {
            content: Box::new(InnerError {
                context,
                ..(*self.content).clone()
            }),
        }
    }

    /// Create a copy of the error with the given underlying errors
    #[must_use]
    pub fn with_underlying_errors(&self, underlying_errors: Vec<Self>) -> Self {
        Self {
            content: Box::new(InnerError {
                underlying_errors,
                ..(*self.content).clone()
            }),
        }
    }

    /// Overwrite the line number with the given number, if applicable
    #[must_use]
    pub fn overwrite_line_number(&self, line_number: usize) -> Self {
        Self {
            content: Box::new(InnerError {
                context: self
                    .content
                    .context
                    .clone()
                    .overwrite_line_number(line_number),
                ..(*self.content).clone()
            }),
        }
    }

    /// Gives the context for this error
    pub const fn context(&self) -> &Context {
        &self.content.context
    }
}

impl fmt::Debug for CustomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "{}: {}{}\n{}",
            self.level(),
            self.content.short_description,
            self.content.context,
            self.content.long_description
        )?;
        match self.content.suggestions.len() {
            0 => Ok(()),
            1 => writeln!(f, "Did you mean: {}?", self.content.suggestions[0]),
            _ => writeln!(
                f,
                "Did you mean any of: {}?",
                self.content.suggestions.join(", ")
            ),
        }?;
        if !self.content.version.is_empty() {
            writeln!(f, "Version: {}", self.content.version)?;
        }
        match self.content.underlying_errors.len() {
            0 => Ok(()),
            1 => writeln!(
                f,
                "Underlying error:\n{}",
                self.content.underlying_errors[0]
            ),
            _ => writeln!(
                f,
                "Underlying errors:\n{}",
                self.content.underlying_errors.iter().join("\n")
            ),
        }
    }
}

impl fmt::Display for CustomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl error::Error for CustomError {}

#[cfg(test)]
#[expect(clippy::print_stdout, clippy::missing_panics_doc)]
mod tests {
    use super::*;
    use crate::error::FilePosition;

    #[test]
    fn create_empty_error() {
        let a = CustomError::error("test", "test", Context::none());
        println!("{a}");
        assert_eq!(format!("{a}"), "error: test\ntest\n");
        assert!(!a.is_warning());
    }

    #[test]
    fn create_full_line_error() {
        let a = CustomError::warning("test", "test", Context::full_line(0, "testing line"));
        println!("{a}");
        assert_eq!(
            format!("{a}"),
            "warning: test\n  ╷\n1 │ testing line\n  ╵\ntest\n"
        );
        assert!(a.is_warning());
    }

    #[test]
    fn create_range_error() {
        let pos1 = FilePosition {
            text: "hello world\nthis is a multiline\npiece of teXt",
            line_index: 0,
            column: 0,
        };
        let pos2 = FilePosition {
            text: "",
            line_index: 3,
            column: 13,
        };
        let a = CustomError::warning("test", "test error", Context::range(&pos1, &pos2));
        println!("{a}");
        assert_eq!(
            format!("{a}"),
            "warning: test\n  ╷\n1 │ hello world\n2 │ this is a multiline\n3 │ piece of teXt\n  ╵\ntest error\n"
        );
        assert!(a.is_warning());
        assert_eq!(pos2.text, "");
        assert_eq!(pos2.line_index, 3);
        assert_eq!(pos2.column, 13);
    }
}
