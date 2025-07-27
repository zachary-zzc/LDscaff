#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
import logging
import random
from pathlib import Path
import vcf
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('LDscaff')

# Path to the C++ matching executable
MATCHING_PROGRAM = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "matching")

def run_command(cmd, description=None, capture_output=True):
    """Run a command and log its output"""
    if description:
        logger.info(description)
    
    logger.debug(f"Command: {' '.join(cmd)}")
    
    if capture_output:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            logger.error(f"Command failed with return code {process.returncode}")
            logger.error(f"STDERR: {stderr}")
            raise RuntimeError(f"Command failed: {' '.join(cmd)}")
        
        logger.debug(f"STDOUT: {stdout}")
        return stdout
    else:
        result = subprocess.run(cmd)
        if result.returncode != 0:
            raise RuntimeError(f"Command failed: {' '.join(cmd)}")

def ensure_directory(directory):
    """Create directory if it doesn't exist"""
    Path(directory).mkdir(parents=True, exist_ok=True)
    return directory

def extract_snp_markers(vcf_file, output_dir, scaffold_names, marker_count=100, minor_ratio=0.0):
    """Extract SNP markers from VCF file for each scaffold"""
    logger.info(f"Extracting SNP markers from {vcf_file}")
    
    # Create output directory
    ensure_directory(output_dir)
    
    # Open VCF reader
    vcf_reader = vcf.Reader(filename=vcf_file)
    
    # Create SNP info log files
    log_f1 = vcf.Writer(open(os.path.join(output_dir, "snps.info.1"), 'w'), vcf_reader)
    log_f2 = vcf.Writer(open(os.path.join(output_dir, "snps.info.2"), 'w'), vcf_reader)
    
    # Sample set
    all_samples = set()
    
    # Extract SNP markers for each scaffold
    for scaffold_name in scaffold_names:
        # Create scaffold directory
        scaffold_dir = os.path.join(output_dir, scaffold_name)
        ensure_directory(scaffold_dir)
        
        logger.info(f"Processing scaffold: {scaffold_name}")
        
        try:
            # Get scaffold length
            scaffold_length = 0
            for record in vcf_reader.fetch(scaffold_name):
                scaffold_length = max(scaffold_length, record.POS)
                
            # Reset reader
            vcf_reader = vcf.Reader(filename=vcf_file)
            
            # Process front end markers (snps.1) - first marker_count markers
            front_snps = {}
            front_records = []
            
            for record in vcf_reader.fetch(scaffold_name, 0, min(5000, scaffold_length // 2)):
                # Check minor allele frequency if needed
                if minor_ratio > 0:
                    valid = check_minor(record, minor_ratio)
                    if not valid:
                        continue
                
                front_records.append(record)
                
                # Extract sample IDs
                for sample in record.samples:
                    all_samples.add(sample.sample)
                    
                if len(front_records) >= marker_count * 2:
                    break
            
            # Take the first marker_count valid records
            for record in front_records[:marker_count]:
                log_f1.write_record(record)
                for sample in record.samples:
                    if sample.sample not in front_snps:
                        front_snps[sample.sample] = [sample.gt_type if sample.gt_type is not None else -1]
                    else:
                        front_snps[sample.sample].append(sample.gt_type if sample.gt_type is not None else -1)
            
            # Write front end markers to files
            for ind, ss in front_snps.items():
                with open(os.path.join(scaffold_dir, f"{ind}.snps.1"), "w") as ofs:
                    ofs.write(" ".join(str(s) for s in ss[:marker_count]))
                    ofs.write("\n")
            
            # Process back end markers (snps.2) - last marker_count markers
            back_snps = {}
            back_records = []
            
            for record in vcf_reader.fetch(scaffold_name, max(0, scaffold_length - 5000), scaffold_length):
                # Check minor allele frequency if needed
                if minor_ratio > 0:
                    valid = check_minor(record, minor_ratio)
                    if not valid:
                        continue
                
                back_records.append(record)
                
                if len(back_records) >= marker_count * 2:
                    break
            
            # Take the last marker_count valid records
            for record in back_records[-marker_count:]:
                log_f2.write_record(record)
                for sample in record.samples:
                    if sample.sample not in back_snps:
                        back_snps[sample.sample] = [sample.gt_type if sample.gt_type is not None else -1]
                    else:
                        back_snps[sample.sample].append(sample.gt_type if sample.gt_type is not None else -1)
            
            # Write back end markers to files
            for ind, ss in back_snps.items():
                with open(os.path.join(scaffold_dir, f"{ind}.snps.2"), "w") as ofs:
                    ofs.write(" ".join(str(s) for s in ss[:marker_count]))
                    ofs.write("\n")
            
        except ValueError as e:
            logger.warning(f"Skipping scaffold {scaffold_name}: {e}")
    
    log_f1.close()
    log_f2.close()
    
    return all_samples

def check_minor(record, cutoff):
    """Check if a record has sufficient minor allele frequency"""
    A = 0
    a = 0
    count = 0
    
    for sample in record.samples:
        gt_type = sample.gt_type
        if gt_type == 0:
            A += 2
            count += 1
        elif gt_type == 1:
            A += 1
            a += 1
            count += 1
        elif gt_type == 2:
            a += 2
            count += 1
        else:
            continue
    
    if count == 0:
        return False
        
    if float(A) / float(2 * count) < cutoff or float(a) / float(2 * count) < cutoff:
        return False
    
    return True

def prepare_samples(output_dir, preprocess_dir, scaffold_names, sample_list):
    """Prepare sample data files for the matching program"""
    logger.info("Preparing sample data files")
    
    # Create samples directory
    samples_dir = os.path.join(output_dir, "samples")
    ensure_directory(samples_dir)
    
    # Prepare files for each sample
    sample_files = []
    for sample in sample_list:
        sample_file = os.path.join(samples_dir, f"{sample}.snps")
        with open(sample_file, 'w') as f_out:
            for scaffold_name in scaffold_names:
                # Process front end markers (snps.1)
                scaffold_dir = os.path.join(preprocess_dir, scaffold_name)
                snp_file_1 = os.path.join(scaffold_dir, f"{sample}.snps.1")
                
                if os.path.exists(snp_file_1):
                    with open(snp_file_1, 'r') as f_in:
                        f_out.write(f_in.read())
                
                # Process back end markers (snps.2)
                snp_file_2 = os.path.join(scaffold_dir, f"{sample}.snps.2")
                if os.path.exists(snp_file_2):
                    with open(snp_file_2, 'r') as f_in:
                        f_out.write(f_in.read())
        
        sample_files.append(sample_file)
    
    # Create sample list file
    sample_list_path = os.path.join(output_dir, "sample.list")
    with open(sample_list_path, 'w') as f:
        for sample_file in sample_files:
            f.write(f"{os.path.abspath(sample_file)}\n")
    
    return sample_list_path, len(sample_files)

def run_matching(sample_list_file, sample_count, scaffold_count, output_file, marker_count=100, minor_ratio=0.2, temp_dir="/tmp", edge_log="edge.log"):
    """Run the C++ matching program"""
    logger.info("Running matching algorithm")
    
    cmd = [
        MATCHING_PROGRAM, 
        "-s", sample_list_file,
        "-n", str(sample_count),
        "-c", str(scaffold_count),
        "-m", str(marker_count),
        "-r", str(minor_ratio),
        "-t", temp_dir,
        "-l", edge_log
    ]
    
    result = run_command(cmd, "Performing maximum weight matching", capture_output=False)
    
    # The matching program produces its output to stdout, which we need to redirect to our output file
    if os.path.exists("matching.out"):
        shutil.copy("matching.out", output_file)
    else:
        logger.warning("Matching output file not found. Check if the program completed successfully.")
    
    return output_file

def find_path(node, edges, threshold):
    """Find a path starting from a node"""
    path = [node]
    
    # Find the matching edge for this node
    matching_edge = None
    for edge in edges:
        if edge[0] == node or edge[1] == node:
            # Skip if this is just connecting the two ends of the same scaffold
            if edge[0] % 2 == 0 and edge[1] - edge[0] == 1:
                continue
            matching_edge = edge
            break
    
    if not matching_edge or matching_edge[2] <= threshold:
        return path
    
    # Add the connected node to the path
    if matching_edge[0] == node:
        next_node = matching_edge[1]
    else:
        next_node = matching_edge[0]
    
    path.append(next_node)
    
    # Add the other end of the scaffold
    if next_node % 2 == 0:
        other_end = next_node + 1
    else:
        other_end = next_node - 1
    
    path.append(other_end)
    
    # Continue extending the path
    while True:
        matching_edge = None
        for edge in edges:
            if (edge[0] == other_end or edge[1] == other_end) and edge[2] > threshold:
                # Skip if this is just connecting the two ends of the same scaffold
                if edge[0] % 2 == 0 and edge[1] - edge[0] == 1:
                    continue
                matching_edge = edge
                break
        
        if not matching_edge:
            break
            
        # Get the next node
        if matching_edge[0] == other_end:
            next_node = matching_edge[1]
        else:
            next_node = matching_edge[0]
            
        # Check if we've already seen this node (cycle detection)
        if next_node in path:
            break
            
        path.append(next_node)
        
        # Add the other end of the scaffold
        if next_node % 2 == 0:
            other_end = next_node + 1
        else:
            other_end = next_node - 1
            
        path.append(other_end)
    
    return path

def find_paths(match_file, scaffold_count, threshold, output_file):
    """Find scaffold paths from matching results"""
    logger.info(f"Finding paths with threshold {threshold}")
    
    # Read edges from matching file
    edges = []
    with open(match_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            edges.append((int(parts[0]), int(parts[1]), float(parts[2])))
    
    # Find all paths
    checked_nodes = []
    paths = []
    
    for i in range(2 * scaffold_count):
        if i not in checked_nodes:
            path = find_path(i, edges, threshold)
            paths.append(path)
            checked_nodes.extend(path)
    
    # Write paths to output file
    with open(output_file, 'w') as f:
        f.write(f"# Paths: {len(paths)}\n")
        f.write(f"# Total edges: {sum([max(0, (len(p)-1)//2) for p in paths])}\n\n")
        
        for i, path in enumerate(paths):
            if len(path) > 1:  # Only write non-trivial paths
                f.write(f"Path {i+1}: {' '.join(map(str, path))}\n")
    
    return paths

def build_contig_names_dict(scaffold_names):
    """Build a dictionary mapping node indices to scaffold names"""
    contig_dict = {}
    for i, name in enumerate(scaffold_names):
        # Map front end
        contig_dict[2*i] = name
        # Map back end
        contig_dict[2*i+1] = name
    
    return contig_dict

def get_contig_name_and_orientation(path, contig_dict):
    """Determine name and orientation of contigs in a path"""
    # Extract scaffold indices (converting node indices to scaffold indices)
    scaffolds = []
    orientations = []
    
    # Process pairs of nodes in the path
    for i in range(0, len(path) - 1, 2):
        # Get the scaffold index for the current node
        scaffold_idx = path[i] // 2
        next_scaffold_idx = path[i+1] // 2
        
        # They should be the same scaffold
        if scaffold_idx != next_scaffold_idx:
            logger.warning(f"Path error: nodes {path[i]} and {path[i+1]} belong to different scaffolds")
            continue
        
        scaffold_name = contig_dict.get(scaffold_idx * 2)
        if not scaffold_name:
            continue
            
        # Determine orientation
        # If path[i] is even (front end) and path[i+1] is odd (back end), orientation is forward
        # If path[i] is odd (back end) and path[i+1] is even (front end), orientation is reverse
        orientation = 1 if path[i] % 2 == 1 else 0  # 1 = reverse, 0 = forward
        
        scaffolds.append(scaffold_name)
        orientations.append(orientation)
    
    return scaffolds, orientations

def assemble_scaffolds(paths, contig_dict, input_fasta, output_fasta, gap_size=100):
    """Assemble scaffolds based on paths"""
    logger.info("Assembling final scaffolds")
    
    # Read input FASTA file
    records = SeqIO.to_dict(SeqIO.parse(input_fasta, "fasta"))
    
    # Assemble scaffolds
    assembled_records = []
    
    for i, path in enumerate(paths):
        if len(path) <= 1:
            continue
            
        # Get scaffold names and orientations
        scaffolds, orientations = get_contig_name_and_orientation(path, contig_dict)
        
        # Skip paths with no scaffolds
        if not scaffolds:
            continue
            
        # Create a new scaffold by joining contigs
        scaffold_seq = ""
        scaffolds_used = []
        
        for scaffold_name, orientation in zip(scaffolds, orientations):
            record = records.get(scaffold_name)
            
            if not record:
                logger.warning(f"Scaffold {scaffold_name} not found in input FASTA")
                continue
                
            # Add gap between scaffolds
            if scaffold_seq:
                scaffold_seq += "N" * gap_size
            
            # Add the sequence in the correct orientation
            if orientation == 1:  # Reverse orientation
                scaffold_seq += str(record.seq.reverse_complement())
            else:  # Forward orientation
                scaffold_seq += str(record.seq)
            
            scaffolds_used.append(scaffold_name)
        
        # Skip if no scaffolds were used
        if not scaffolds_used:
            continue
        
        # Create a new record for the scaffold
        scaffold_record = SeqRecord(
            Seq(scaffold_seq),
            id=f"scaffold_{i+1}",
            description=f"Assembled from {', '.join(scaffolds_used)}"
        )
        
        assembled_records.append(scaffold_record)
    
    # Write output FASTA file
    SeqIO.write(assembled_records, output_fasta, "fasta")
    logger.info(f"Wrote {len(assembled_records)} scaffolds to {output_fasta}")

def get_scaffold_names_from_fasta(fasta_file):
    """Extract scaffold names from FASTA file"""
    scaffold_names = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Extract scaffold name (removing '>' and any description after whitespace)
                scaffold_name = line[1:].split()[0].strip()
                scaffold_names.append(scaffold_name)
    
    return scaffold_names

def get_scaffold_names_from_vcf(vcf_file):
    """Extract scaffold names from VCF file"""
    vcf_reader = vcf.Reader(filename=vcf_file)
    scaffold_names = set()
    
    # Try to get contigs from header
    if hasattr(vcf_reader, 'contigs'):
        for contig in vcf_reader.contigs.values():
            scaffold_names.add(contig.id)
    
    # If no contigs in header, get from records
    if not scaffold_names:
        for record in vcf_reader:
            scaffold_names.add(record.CHROM)
            # Limit to 100 records to avoid scanning the whole file
            if len(scaffold_names) >= 100:
                break
    
    return list(scaffold_names)

def run_full_pipeline(args):
    """Run the complete LDscaff pipeline"""
    # Create main output directory
    output_dir = ensure_directory(args.output)
    
    # Get scaffold names from FASTA file
    if args.input_fasta:
        scaffold_names = get_scaffold_names_from_fasta(args.input_fasta)
    else:
        # Get scaffold names from VCF file
        scaffold_names = get_scaffold_names_from_vcf(args.vcf)
    
    logger.info(f"Found {len(scaffold_names)} scaffolds")
    
    # Extract SNP markers
    preprocess_dir = os.path.join(output_dir, "preprocess")
    samples = extract_snp_markers(
        args.vcf,
        preprocess_dir,
        scaffold_names,
        args.marker_count,
        args.minor_ratio
    )
    
    # Get or create sample list
    if args.samples:
        with open(args.samples, 'r') as f:
            sample_list = [line.strip() for line in f if line.strip()]
    else:
        sample_list = list(samples)
        sample_list_file = os.path.join(output_dir, "samples.txt")
        with open(sample_list_file, 'w') as f:
            for sample in sample_list:
                f.write(f"{sample}\n")
    
    # Prepare sample data for matching
    sample_list_file, sample_count = prepare_samples(output_dir, preprocess_dir, scaffold_names, sample_list)
    
    # Run matching
    edge_log = os.path.join(output_dir, "edge.log")
    match_result = os.path.join(output_dir, "match.result")
    
    run_matching(
        sample_list_file, 
        sample_count, 
        len(scaffold_names), 
        match_result,
        args.marker_count,
        args.minor_ratio,
        args.temp_dir,
        edge_log
    )
    
    # Find paths
    paths_result = os.path.join(output_dir, "paths.result")
    paths = find_paths(match_result, len(scaffold_names), args.ld_threshold, paths_result)
    
    # Build scaffold dictionary
    scaffold_dict = {}
    for i, name in enumerate(scaffold_names):
        scaffold_dict[i] = name
    
    # Assemble scaffolds if input FASTA provided
    if args.input_fasta:
        output_fasta = args.output_fasta if args.output_fasta else os.path.join(output_dir, "assembled.fasta")
        assemble_scaffolds(paths, scaffold_dict, args.input_fasta, output_fasta, args.gap_size)
    
    logger.info("LDscaff pipeline completed successfully")

def cmd_extract_markers(args):
    """Command to extract SNP markers from VCF file"""
    # Get scaffold names
    if args.scaffolds:
        with open(args.scaffolds, 'r') as f:
            scaffold_names = [line.strip() for line in f if line.strip()]
    else:
        scaffold_names = get_scaffold_names_from_vcf(args.vcf)
    
    ensure_directory(args.output)
    samples = extract_snp_markers(
        args.vcf,
        args.output,
        scaffold_names,
        args.marker_count,
        args.minor_ratio
    )
    
    # Write sample list if requested
    if args.sample_list:
        with open(args.sample_list, 'w') as f:
            for sample in samples:
                f.write(f"{sample}\n")
        
    logger.info(f"Extracted markers for {len(scaffold_names)} scaffolds and {len(samples)} samples")

def cmd_prepare_samples(args):
    """Command to prepare sample data for matching"""
    ensure_directory(args.output)
    
    # Get scaffold names
    if args.scaffolds:
        with open(args.scaffolds, 'r') as f:
            scaffold_names = [line.strip() for line in f if line.strip()]
    else:
        # Get scaffold names from marker directories
        scaffold_names = []
        for item in os.listdir(args.marker_dir):
            if os.path.isdir(os.path.join(args.marker_dir, item)):
                scaffold_names.append(item)
    
    # Read sample list
    with open(args.samples, 'r') as f:
        sample_list = [line.strip() for line in f if line.strip()]
    
    # Prepare samples
    sample_list_file, sample_count = prepare_samples(args.output, args.marker_dir, scaffold_names, sample_list)
    logger.info(f"Prepared sample data files for {sample_count} samples")
    logger.info(f"Sample list file: {sample_list_file}")
    
    # Write scaffold and sample counts
    with open(os.path.join(args.output, "counts.txt"), 'w') as f:
        f.write(f"Scaffolds: {len(scaffold_names)}\n")
        f.write(f"Samples: {sample_count}\n")

def cmd_matching(args):
    """Command to run matching algorithm"""
    # Ensure output directory exists
    ensure_directory(os.path.dirname(args.output))
    
    # Run matching
    run_matching(
        args.sample_list, 
        args.samples, 
        args.scaffolds, 
        args.output,
        args.marker_count,
        args.minor_ratio,
        args.temp_dir,
        args.edge_log
    )
    logger.info(f"Matching completed. Results written to {args.output}")

def cmd_find_paths(args):
    """Command to find scaffold paths"""
    # Ensure output directory exists
    ensure_directory(os.path.dirname(args.output))
    
    # Find paths
    paths = find_paths(args.match_result, args.scaffolds, args.threshold, args.output)
    logger.info(f"Found {len(paths)} paths")
    logger.info(f"Wrote path results to {args.output}")

def cmd_assemble(args):
    """Command to assemble scaffolds"""
    # Ensure output directory exists
    ensure_directory(os.path.dirname(args.output_fasta))
    
    # Read paths
    paths = []
    with open(args.paths, 'r') as f:
        for line in f:
            if line.startswith("Path"):
                parts = line.strip().split(": ")
                if len(parts) > 1:
                    path = [int(x) for x in parts[1].split()]
                    paths.append(path)
    
    # Get scaffold names
    if args.scaffolds:
        with open(args.scaffolds, 'r') as f:
            scaffold_names = [line.strip() for line in f if line.strip()]
    else:
        scaffold_names = get_scaffold_names_from_fasta(args.input_fasta)
    
    # Build scaffold dictionary
    scaffold_dict = {}
    for i, name in enumerate(scaffold_names):
        scaffold_dict[i] = name
    
    # Assemble scaffolds
    assemble_scaffolds(paths, scaffold_dict, args.input_fasta, args.output_fasta, args.gap_size)

def main():
    """Main entry point for the CLI"""
    parser = argparse.ArgumentParser(
        description="LDscaff: LD-based scaffolding of de novo genome assemblies",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Global arguments
    parser.add_argument('--verbose', action='store_true', 
                        help='Enable verbose output')
    
    # Create subparsers for each command
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Full pipeline command
    pipeline_parser = subparsers.add_parser('run', help='Run the complete LDscaff pipeline')
    pipeline_parser.add_argument('--vcf', required=True,
                        help='Input VCF file with population variation data')
    pipeline_parser.add_argument('--output', required=True, 
                        help='Output directory')
    pipeline_parser.add_argument('--samples', 
                        help='Sample list file (optional, will be generated if not provided)')
    pipeline_parser.add_argument('--gap-size', type=int, default=100, 
                        help='Gap size between scaffolds in final assembly')
    pipeline_parser.add_argument('--minor-ratio', type=float, default=0.2, 
                        help='SNP marker minor ratio cutoff')
    pipeline_parser.add_argument('--ld-threshold', type=float, default=0.2, 
                        help='LD threshold for breaking weak connections')
    pipeline_parser.add_argument('--marker-count', type=int, default=100,
                        help='Number of markers per scaffold end')
    pipeline_parser.add_argument('--input-fasta', 
                        help='Input scaffold FASTA file for final assembly (optional)')
    pipeline_parser.add_argument('--output-fasta', 
                        help='Output assembled scaffold FASTA file (optional)')
    pipeline_parser.add_argument('--temp-dir', default="/tmp",
                        help='Temporary directory for matching program')
    pipeline_parser.set_defaults(func=run_full_pipeline)
    
    # Extract markers command
    markers_parser = subparsers.add_parser('extract', help='Extract SNP markers from VCF file')
    markers_parser.add_argument('--vcf', required=True, 
                        help='Input VCF file')
    markers_parser.add_argument('--output', required=True, 
                        help='Output directory for markers')
    markers_parser.add_argument('--scaffolds', 
                        help='File containing scaffold names (optional, will be extracted from VCF if not provided)')
    markers_parser.add_argument('--marker-count', type=int, default=100, 
                        help='Number of markers per scaffold end')
    markers_parser.add_argument('--minor-ratio', type=float, default=0.2, 
                        help='SNP marker minor ratio cutoff')
    markers_parser.add_argument('--sample-list', 
                        help='Output file for sample list (optional)')
    markers_parser.set_defaults(func=cmd_extract_markers)
    
    # Prepare samples command
    prepare_parser = subparsers.add_parser('prepare', help='Prepare sample data for matching')
    prepare_parser.add_argument('--marker-dir', required=True, 
                        help='Directory containing extracted markers')
    prepare_parser.add_argument('--scaffolds', 
                        help='File containing scaffold names (optional, will be extracted from marker directory if not provided)')
    prepare_parser.add_argument('--samples', required=True, 
                        help='Sample list file')
    prepare_parser.add_argument('--output', required=True, 
                        help='Output directory')
    prepare_parser.set_defaults(func=cmd_prepare_samples)
    
    # Matching command
    matching_parser = subparsers.add_parser('match', help='Run matching algorithm')
    matching_parser.add_argument('--sample-list', required=True, 
                        help='Sample list file from prepare step')
    matching_parser.add_argument('--samples', required=True, type=int, 
                        help='Number of samples')
    matching_parser.add_argument('--scaffolds', required=True, type=int, 
                        help='Number of scaffolds')
    matching_parser.add_argument('--output', required=True, 
                        help='Output file for matching results')
    matching_parser.add_argument('--marker-count', type=int, default=100,
                        help='Number of markers to use from each scaffold')
    matching_parser.add_argument('--minor-ratio', type=float, default=0.2,
                        help='Minor allele cutoff')
    matching_parser.add_argument('--temp-dir', default="/tmp",
                        help='Temporary directory')
    matching_parser.add_argument('--edge-log', default="edge.log",
                        help='Log file for edge weights')
    matching_parser.set_defaults(func=cmd_matching)
    
    # Find paths command
    paths_parser = subparsers.add_parser('paths', help='Find scaffold paths')
    paths_parser.add_argument('--match-result', required=True, 
                        help='Match result file from match step')
    paths_parser.add_argument('--scaffolds', required=True, type=int, 
                        help='Number of scaffolds')
    paths_parser.add_argument('--threshold', required=True, type=float, 
                        help='LD threshold for breaking weak connections')
    paths_parser.add_argument('--output', required=True, 
                        help='Output file for path results')
    paths_parser.set_defaults(func=cmd_find_paths)
    
    # Assemble command
    assemble_parser = subparsers.add_parser('assemble', help='Assemble scaffolds')
    assemble_parser.add_argument('--paths', required=True, 
                        help='Path result file from paths step')
    assemble_parser.add_argument('--scaffolds', 
                        help='File containing scaffold names (optional, will be extracted from input FASTA if not provided)')
    assemble_parser.add_argument('--input-fasta', required=True, 
                        help='Input scaffold FASTA file')
    assemble_parser.add_argument('--output-fasta', required=True, 
                        help='Output assembled scaffold FASTA file')
    assemble_parser.add_argument('--gap-size', type=int, default=100, 
                        help='Gap size between scaffolds in output')
    assemble_parser.set_defaults(func=cmd_assemble)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Run the appropriate command
    if args.command is None:
        parser.print_help()
    else:
        try:
            args.func(args)
        except Exception as e:
            logger.error(f"Error in LDscaff: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)

if __name__ == "__main__":
    main()
