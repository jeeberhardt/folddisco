use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::Path;

use flate2::read::GzDecoder;

use crate::structure::atom::Atom;

use super::super::core::*;
use super::*;

/// A PDB reader
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    /// The underlying reader
    pub reader: R,
    ///
    pub input_type: StructureFileFormat,
}

// Add a new error type for PDB parsing
#[derive(Debug)]
pub struct PDBParseError {
    pub message: String,
}

impl std::fmt::Display for PDBParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for PDBParseError {}

// Helper to extract residue name as 3-char array, accepting numbers as valid
fn parse_residue_name(resname: &str) -> [u8; 3] {
    let bytes = resname.as_bytes();
    match bytes.len() {
        1 => [bytes[0], b' ', b' '],
        2 => [bytes[0], bytes[1], b' '],
        3 => [bytes[0], bytes[1], bytes[2]],
        _ => [b' ', b' ', b' '],
    }
}

// Helper to extract atom name as 4-char array, accepting short names as right-justified
fn parse_atom_name(atom_name: &str) -> [u8; 4] {
    let bytes = atom_name.as_bytes();
    match bytes.len() {
        1 => [b' ', b' ', b' ', bytes[0]],
        2 => [b' ', b' ', bytes[0], bytes[1]],
        3 => [b' ', bytes[0], bytes[1], bytes[2]],
        4 => [bytes[0], bytes[1], bytes[2], bytes[3]],
        _ => [b' ', b' ', b' ', b' '],
    }
}

// Update parse_line to propagate errors and use consistent error messages
fn parse_line(line: &str) -> Result<Atom, PDBParseError> {
    // Example PDB ATOM line format:
    // COLUMNS        DATA TYPE       FIELD         DEFINITION
    // -------------------------------------------------------------------------------------
    // 1 - 6         Record name     "ATOM  "
    // 7 - 11        Integer         serial        Atom serial number.
    // 13 - 16       Atom            name          Atom name.
    // 17            Character       altLoc        Alternate location indicator.
    // 18 - 20       Residue name    resName       Residue name.
    // 22            Character       chainID       Chain identifier.
    // 23 - 26       Integer         resSeq        Residue sequence number.
    // 31 - 38       Real(8.3)       x             Orthogonal coordinates for X in Angstroms.
    // 39 - 46       Real(8.3)       y             Orthogonal coordinates for Y in Angstroms.
    // 47 - 54       Real(8.3)       z             Orthogonal coordinates for Z in Angstroms.
    // 55 - 60       Real(6.2)       occupancy     Occupancy.
    // 61 - 66       Real(6.2)       tempFactor    Temperature factor.
    // 77 - 78       LString(2)      element       Element symbol, right-justified.
    // 79 - 80       LString(2)      charge        Charge on the atom.

    let atom_name = line.get(12..16).unwrap_or("").trim();
    let residue_name = line.get(17..20).unwrap_or("").trim();
    let chain_id = line.get(21..22).unwrap_or("").trim();
    let residue_number = line.get(22..26).unwrap_or("").trim();
    let x = line.get(30..38).unwrap_or("").trim();
    let y = line.get(38..46).unwrap_or("").trim();
    let z = line.get(46..54).unwrap_or("").trim();
    let b_factor = line.get(60..66).unwrap_or("").trim();

    // Consistent error reporting
    let attempted_atom_name = if !atom_name.is_empty() { atom_name } else { "missing" };
    let attempted_residue_name = if !residue_name.is_empty() { residue_name } else { "missing" };

    if atom_name.is_empty() {
        return Err(PDBParseError {
            message: format!(
                "Atom name missing in atom record (atom: '{}', residue: '{}')",
                attempted_atom_name, attempted_residue_name
            ),
        });
    }
    if residue_name.is_empty() {
        return Err(PDBParseError {
            message: format!(
                "Residue name missing in atom record (atom: '{}', residue: '{}')",
                attempted_atom_name, attempted_residue_name
            ),
        });
    }
    if residue_number.is_empty() {
        return Err(PDBParseError {
            message: format!(
                "Residue number missing in atom record (atom: '{}', residue: '{}')",
                attempted_atom_name, attempted_residue_name
            ),
        });
    }
    if chain_id.is_empty() {
        return Err(PDBParseError {
            message: format!(
                "Chain name missing in atom record (atom: '{}', residue: '{}')",
                attempted_atom_name, attempted_residue_name
            ),
        });
    }
    if x.is_empty() {
        return Err(PDBParseError {
            message: format!(
                "Atom X position missing in atom record (atom: '{}', residue: '{}')",
                attempted_atom_name, attempted_residue_name
            ),
        });
    }
    if y.is_empty() {
        return Err(PDBParseError {
            message: format!(
                "Atom Y position missing in atom record (atom: '{}', residue: '{}')",
                attempted_atom_name, attempted_residue_name
            ),
        });
    }
    if z.is_empty() {
        return Err(PDBParseError {
            message: format!(
                "Atom Z position missing in atom record (atom: '{}', residue: '{}')",
                attempted_atom_name, attempted_residue_name
            ),
        });
    }

    // Parse numeric fields
    let id = line.get(6..11).unwrap_or("").trim().parse::<u64>().unwrap_or(0);
    let residue_number = residue_number.parse::<u64>().unwrap_or(0);
    let pos_x = x.parse::<f32>().unwrap_or(0.0);
    let pos_y = y.parse::<f32>().unwrap_or(0.0);
    let pos_z = z.parse::<f32>().unwrap_or(0.0);
    let b_factor = b_factor.parse::<f32>().unwrap_or(1.0);
    let residue_name_arr = parse_residue_name(residue_name);
    let chain_name_vec = chain_id.as_bytes().to_vec();

    Ok(Atom::new(
        pos_x, pos_y, pos_z,
        parse_atom_name(atom_name),
        id,
        chain_name_vec,
        residue_name_arr,
        residue_number,
        b_factor,
    ))
}

// Update read_structure and read_structure_from_gz to propagate and print errors
impl Reader<File> {
    pub fn new(file: File) -> Self {
        Reader {
            reader: file,
            input_type: StructureFileFormat::PDB,
        }
    }

    /// Read from a file path
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> Result<Self, &'static str> {
        File::open(&path)
            .map(Reader::new)
            .map_err(|_e| "Error opening file")
    }

    pub fn read_structure(&self) -> Result<Structure, &str> {
        let reader = BufReader::new(&self.reader);
        let mut structure = Structure::new(); // revise
        let mut record = (vec![b' '], 0);
        let mut model = 0;
        let mut errors = Vec::new();
        // Reading each line of PDB, parse and build atomvector.
        for (_idx, line) in reader.lines().enumerate() {
            if let Ok(atomline) = line {
                if model > 1 {
                    // Current version does not support multiple models in one PDB file
                    break;
                }
                // If line is less than 6 characters, skip the line
                if atomline.len() < 6 {
                    continue;
                }
                match &atomline[..6] {
                    "MODEL " => {
                        model += 1;
                    }
                    "ATOM  " => {
                        match parse_line(&atomline) {
                            Ok(atom) => {
                                structure.update(atom, &mut record);
                            }
                            Err(e) => {
                                errors.push(e);
                                continue;
                            }
                        }
                    }
                    _ => continue,
                }
            } else {
                return Err("Error reading line");
            };
        }
        if !errors.is_empty() {
            let file_name = "<unknown>"; // You can add a path field to Reader if you want
            eprintln!("\nError(s) in file: {}\n", file_name);
            for error in &errors {
                eprintln!("{}", error);
            }
            return Err("Error parsing PDB file");
        }
        // println!("{structure:?}");
        Ok(structure)
    }

    pub fn read_structure_from_gz(&self) -> Result<Structure, &str> {
        // Load whole file and close the file
        let mut decoder = GzDecoder::new(&self.reader);
        let mut binary = Vec::new();
        decoder.read_to_end(&mut binary).unwrap();
        // Close the decoder after flushing
        decoder.flush().unwrap();
        // Drop the decoder
        drop(decoder);

        // Create a new Structure
        let mut structure = Structure::new();
        let mut record = (vec![b' '], 0);
        let reader = BufReader::new(&binary[..]);
        let mut errors = Vec::new();
        
        // Read binary as a string. Conver
        for (_idx, line) in reader.lines().enumerate() {
            if let Ok(atomline) = line {
                match &atomline[..6] {
                    "ATOM  " => {
                        match parse_line(&atomline) {
                            Ok(atom) => {
                                structure.update(atom, &mut record);
                            }
                            Err(e) => {
                                errors.push(e);
                                continue;
                            }
                        }
                    }
                    _ => continue,
                }
            } else {
                return Err("Error reading line");
            };
        }
        if !errors.is_empty() {
            let file_name = "<unknown>"; // You can add a path field to Reader if you want
            eprintln!("\nError(s) in file: {}\n", file_name);
            for error in &errors {
                eprintln!("{}", error);
            }
            return Err("Error parsing PDB file");
        }
        // Drop the binary
        drop(binary);
        Ok(structure)
    }
    
}

#[cfg(test)]
mod tests {
    use crate::prelude::load_path;

    use super::*;
    use std::fs::File;
    
    use std::path::Path;
    
    
    #[test]
    fn test_read_pdb() {
        let path = Path::new("data/homeobox/1akha-.pdb");
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure().unwrap();
        let compact = structure.to_compact();
        assert_eq!(compact.num_residues, 49);
    }
    #[test]
    fn test_get_min_max_coords() {
        let path = Path::new("data/homeobox/1akha-.pdb");
        // let path = Path::new("analysis/h_sapiens_pdb/AF-Q02817-F3-model_v4.pdb");
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure().unwrap();
        let compact = structure.to_compact();
        let min = compact.ca_vector.min_coord();
        let max = compact.ca_vector.max_coord();
        let min = min.unwrap();
        let max = max.unwrap();
        let len = compact.ca_vector.x.len();
        // Distance between min and max
        let dist = max.distance(&min);
        let diff = max.sub(&min);
        let mut bins = diff.clone();
        bins.x = (bins.x / 30.0).ceil().min(6.0);
        bins.y = (bins.y / 30.0).ceil().min(6.0);
        bins.z = (bins.z / 30.0).ceil().min(6.0);
        println!("Bins: {:?}", bins);
        // let dist = (dist.0.powi(2) + dist.1.powi(2) + dist.2.powi(2)).sqrt();
        let x_bin = diff.x / bins.x;
        let y_bin = diff.y / bins.y;
        let z_bin = diff.z / bins.z;
        println!("Min: {:?}, Max: {:?}, len: {}, dist: {:?}, x_bin: {}, y_bin: {}, z_bin: {}", min, max, len, dist, x_bin, y_bin, z_bin);
    }
    
    
    
    #[test]
    fn test_read_pdb_gz() {
        let path = Path::new("data/homeobox/inner/1akha-.pdb.gz");
        let file = File::open(&path).unwrap();
        let reader = Reader::new(file);
        let structure = reader.read_structure_from_gz().unwrap();
        let compact = structure.to_compact();
        assert_eq!(compact.num_residues, 49);
    }
    
    #[test]
    fn test_loading_works() {
        let dir = "data/io_test";
        let pdb_paths = load_path(dir, false);
        for pdb_path in pdb_paths {
            let file = File::open(&pdb_path).unwrap();
            let reader = Reader::new(file);
            // If the file is gzipped, use the gzipped reader
            let structure = if pdb_path.ends_with(".gz") {
                reader.read_structure_from_gz().unwrap()
            } else {
                reader.read_structure().unwrap()
            };
            let compact = structure.to_compact();
            println!("{}:{}", pdb_path, compact.num_residues);
        }
    }
}