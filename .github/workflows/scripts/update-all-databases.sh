#!/usr/bin/env bash


function help {
    echo "Usage: generate-all-databases.sh"
    echo ""
    echo "Download the required databases and build the "
    echo "required binary representations of the ontologies."
    echo ""
    echo "Options:"
    echo "  -h, --help     Display this help and exit"
    exit 1
}


# Download IMGT and process and serialize it to a binary blob.
function make-imgt {
    echo "Downloading IMGT..."
    mkdir -p rustyms-generate-imgt/data
    # IMGT is not very reliable, so sometimes the server is down.
    # The || clause here allows the rest of the script to continue
    # even if this fails.
    curl https://www.imgt.org/download/LIGM-DB/imgt.dat.Z \
        | gunzip -c > rustyms-generate-imgt/data/imgt.dat \
        && echo "Serializing IMGT ..." \
        && cargo run -p rustyms-generate-imgt \
        || echo "Failed to download IMGT. I did not update it." >> /tmp/MESSAGES
}


# Download the relevant ontologies and serialize them to binary blobs.
function make-ontologies {
    echo "Downloading databases..."
    db_data="rustyms-generate-databases/data"
    mkdir -p ${db_data}
    curl https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/refs/heads/master/PSI-MOD-newstyle.obo \
        > ${db_data}/PSI-MOD-newstyle.obo
    curl https://unimod.org/xml/unimod.xml > ${db_data}/unimod.xml
    curl ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/RESIDUES.XML \
        > ${db_data}/RESID-RESIDUES.XML
    curl https://raw.githubusercontent.com/HUPO-PSI/mzIdentML/master/cv/XLMOD.obo \
        > ${db_data}/XLMOD.obo
    curl -L http://purl.obolibrary.org/obo/gno.obo \
        > ${db_data}/GNOme.obo
    curl -L https://glycosmos.org/download/glycosmos_glycans_list.csv \
        > ${db_data}/glycosmos_glycans_list.csv


    echo "Serializing the other databases..."
    cargo run -p rustyms-generate-databases
}


function main {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -h|--help)
                help
                ;;
            *)
                echo "Unknown argument: $1"
                help
                ;;
        esac
    done

    touch /tmp/MESSAGES

    make-imgt
    make-ontologies
}


main "$@"
