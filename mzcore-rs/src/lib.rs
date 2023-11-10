pub mod chemistry;
pub mod common;
pub mod ms;
pub mod msms;

#[cfg(test)]
mod tests {
    use crate::ms::utils::MassTolWindow;

    #[test]
    fn tolerances() {
        assert_eq!(
            MassTolWindow::ppm(-10.0, 20.0).bounds(1000.0),
            (999.99, 1000.02)
        );
    }
}
