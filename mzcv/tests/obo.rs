//! Parse all Obo files to test the parser

use mzcv::OboOntology;

#[ignore = "Only run when interested in extending the support of the Obo parser (example files are not included)"]
#[test]
fn test_obo_files() {
    for file in std::fs::read_dir("data").unwrap() {
        let file = file.unwrap();
        if file
            .path()
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("obo"))
        {
            let _obo = OboOntology::from_file(file.path()).unwrap();
        }
    }
}
