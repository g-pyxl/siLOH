import sys
import json
import csv
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional, Set
from pathlib import Path
import logging
from enum import Enum

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

class Sex(Enum):
    MALE = 'Male'
    FEMALE = 'Female'
    UNKNOWN = 'Unknown'

@dataclass
class GenomicRegion:
    chromosome: str
    start: int
    end: int
    homozygous_count: int
    total_count: int
    
    @property
    def size(self) -> int:
        return self.end - self.start + 1
    
    @property
    def confidence(self) -> float:
        return (self.homozygous_count / self.total_count * 100) if self.total_count > 0 else 0

class LOHAnalyzer:
    def __init__(
        self,
        min_streak: int = 5,
        loh_threshold: float = 35.0,
        min_region_size: int = 1_000_000,
        max_gap: int = 2,
        sex_determination_threshold: float = 0.2
    ):
        self.min_streak = min_streak
        self.loh_threshold = loh_threshold
        self.min_region_size = min_region_size
        self.max_gap = max_gap
        self.sex_determination_threshold = sex_determination_threshold
        
    def _process_position(
        self, 
        chrom: str, 
        pos: int, 
        var_freq: float,
        current_region: Optional[GenomicRegion],
        gap_count: int,
        regions: List[GenomicRegion]
    ) -> Tuple[Optional[GenomicRegion], int]:
        """Process a single genomic position and update LOH regions."""
        is_homozygous = self._is_homozygous(var_freq)
        
        # Handle chromosome change
        if current_region and current_region.chromosome != chrom:
            if current_region:
                regions.append(current_region)
            current_region = None
            gap_count = 0
        
        # Process homozygous position
        if is_homozygous:
            if current_region is None:
                current_region = GenomicRegion(chrom, pos, pos, 1, 1)
            else:
                current_region = GenomicRegion(
                    chrom,
                    current_region.start,
                    pos,
                    current_region.homozygous_count + 1,
                    current_region.total_count + 1
                )
            gap_count = 0
        
        # Process heterozygous position
        else:
            if current_region is not None:
                gap_count += 1
                current_region = GenomicRegion(
                    chrom,
                    current_region.start,
                    pos,
                    current_region.homozygous_count,
                    current_region.total_count + 1
                )
                
                # Close region if gap is too large
                if gap_count > self.max_gap:
                    regions.append(current_region)
                    current_region = None
                    gap_count = 0
        
        return current_region, gap_count

    @staticmethod
    def load_centromeres(json_path: Path) -> Dict[str, int]:
        """Load centromere positions from JSON file."""
        try:
            with open(json_path) as f:
                data = json.load(f)
            return {k: v['centromere'] for k, v in data.items()}
        except Exception as e:
            logging.error(f"Failed to load centromeres file: {e}")
            raise

    @staticmethod
    def load_bed_regions(bed_path: Path) -> Dict[str, List[Tuple[int, int, str]]]:
        """Load gene regions from BED file."""
        regions = {}
        try:
            with open(bed_path) as f:
                for line in f:
                    chrom, start, end, gene = line.strip().split('\t')
                    chrom = f"chr{chrom}" if not chrom.startswith('chr') else chrom
                    regions.setdefault(chrom, []).append((int(start), int(end), gene))
            return regions
        except Exception as e:
            logging.error(f"Failed to load BED file: {e}")
            raise

    def _is_homozygous(self, var_freq: float) -> bool:
        """Determine if a variant frequency indicates homozygosity."""
        return var_freq <= self.loh_threshold or var_freq >= (100 - self.loh_threshold)

    def analyze_file(self, file_path: Path, centromeres: Dict[str, int]) -> Tuple[List[GenomicRegion], Sex]:
        """Analyze a CNS file for LOH regions and determine sample sex."""
        regions: List[GenomicRegion] = []
        current_region: Optional[GenomicRegion] = None
        gap_count = 0
        chrX_stats = {'het': 0, 'total': 0}
        
        try:
            with open(file_path) as f:
                next(f)  # Skip header
                for line in f:
                    chrom, pos, *_, var_freq = line.strip().split('\t')[:7]
                    pos = int(pos)
                    var_freq = float(var_freq.strip('%'))
                    
                    # Track chrX heterozygosity
                    if chrom == 'chrX':
                        chrX_stats['total'] += 1
                        if not self._is_homozygous(var_freq):
                            chrX_stats['het'] += 1
                    
                    # Process LOH regions
                    current_region, gap_count = self._process_position(
                        chrom, pos, var_freq, current_region, gap_count, regions
                    )
            
            # Add the last region if it exists
            if current_region:
                regions.append(current_region)
            
            # Determine sex
            sex = self._determine_sex(chrX_stats['het'], chrX_stats['total'])
            
            return self._filter_regions(regions, centromeres), sex
            
        except Exception as e:
            logging.error(f"Error processing file {file_path}: {e}")
            raise

    def _determine_sex(self, het_count: int, total_count: int) -> Sex:
        """Determine sample sex based on X chromosome heterozygosity."""
        if total_count == 0:
            return Sex.UNKNOWN
        ratio = het_count / total_count
        return Sex.FEMALE if ratio > self.sex_determination_threshold else Sex.MALE

    def _filter_regions(self, regions: List[GenomicRegion], centromeres: Dict[str, int]) -> List[GenomicRegion]:
        """Filter and split regions based on size and centromere positions."""
        filtered_regions = []
        for region in regions:
            if (region.homozygous_count >= self.min_streak and 
                region.size >= self.min_region_size):
                centromere_pos = centromeres.get(region.chromosome)
                if centromere_pos and region.start < centromere_pos < region.end:
                    # Split region at centromere
                    filtered_regions.extend([
                        GenomicRegion(
                            region.chromosome, region.start, centromere_pos - 1,
                            region.homozygous_count // 2, region.total_count // 2
                        ),
                        GenomicRegion(
                            region.chromosome, centromere_pos, region.end,
                            region.homozygous_count // 2, region.total_count // 2
                        )
                    ])
                else:
                    filtered_regions.append(region)
        return filtered_regions

    def find_affected_genes(self, region: GenomicRegion, bed_regions: Dict[str, List[Tuple[int, int, str]]]) -> Set[str]:
        """Find unique genes that overlap with a given genomic region."""
        affected_genes = set()  # Changed from list to set
        for start, end, gene in bed_regions.get(region.chromosome, []):
            if start <= region.end and end >= region.start:
                affected_genes.add(gene)  # Using add instead of append
        return affected_genes

class ResultsWriter:
    @staticmethod
    def save_to_csv(output_path: Path, results: List[Tuple[str, int, int, str]]):
        """Save analysis results to CSV file."""
        try:
            with open(output_path, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, 
                    fieldnames=['Chromosome', 'Start', 'End', 'Affected_Genes'])
                writer.writeheader()
                for chrom, start, end, genes in results:
                    writer.writerow({
                        'Chromosome': chrom,
                        'Start': start,
                        'End': end,
                        'Affected_Genes': genes
                    })
        except Exception as e:
            logging.error(f"Failed to save results to {output_path}: {e}")
            raise

def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        logging.error("Usage: python floh.py <file_path> [<bed_file>]")
        sys.exit(1)

    try:
        analyzer = LOHAnalyzer()
        file_path = Path(sys.argv[1])
        bed_path = Path(sys.argv[2]) if len(sys.argv) == 3 else None
        centromeres_path = Path('centromeres.json')

        # Load required files
        centromeres = analyzer.load_centromeres(centromeres_path)
        bed_regions = analyzer.load_bed_regions(bed_path) if bed_path else None

        # Analyze file
        loh_regions, sex = analyzer.analyze_file(file_path, centromeres)

        # Process results
        results = []
        for region in loh_regions:
            if 'X' not in region.chromosome and region.homozygous_count > 40 and region.confidence > 90:
                if bed_regions:
                    affected_genes = analyzer.find_affected_genes(region, bed_regions)
                    if affected_genes:  # Only include regions with affected genes
                        results.append((
                            region.chromosome,
                            region.start,
                            region.end,
                            ','.join(sorted(affected_genes))  # Sort genes for consistent output
                        ))
                else:
                    results.append((
                        region.chromosome,
                        region.start,
                        region.end,
                        ''  # No genes information when BED file is not provided
                    ))

        # Save results
        output_path = file_path.with_suffix('.loh.csv')
        ResultsWriter.save_to_csv(output_path, results)

        # Log results
        logging.info(f"Analysis complete. Sex: {sex.value}")
        logging.info(f"Results saved to {output_path}")

    except Exception as e:
        logging.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()