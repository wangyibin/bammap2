use std::collections::HashMap;
use calign::index::index;


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index() {
        
        let ref_genome: Vec<_> = vec![("seq1".to_string(), b"CGGCAmmAGGTTAAAATCTmAGTGCTGCAATAGGCGATTACAGTACAGCACCCAGCCTCCC".to_vec())].into_iter().collect();
        // ref_genome to HashMap 
        let ref_genome = ref_genome.iter().map(|(k, v)| (k.clone(), v.clone())).collect::<HashMap<String, Vec<u8>>>();
        let k = 15;
        let w = 5;
        let index = index(&ref_genome, k, w);
    
    }
}