import gzip

# We use this for reverse complementing and we don't 
# want to remake it every time.
revcomp_table = str.maketrans("GATC","CTAG")
stop_codons = {"TAA","TAG","TGA"}

class CodingTransposons:

    def __init__ (self, transposons, hits):
        self.transposons = transposons
        self.hits = hits


    def parseHits(self):

        # We'll save the hits by chromosome
        self.coding_transposons = {}

        with gzip.open(self.hits, "rt", encoding="utf8") as infh:
            for line in infh:
                if line.startswith("#"):
                    continue
                sections = line.split("\t")
                if not sections[1] in self.transposons:
                    # This isn't a transposon which we have annotated with a CDS
                    continue
                
                # We need to find if this covers a CDS

                # The positions within the repeat consensus are always
                # on the positive strand, so start is lower than end
                hit_start = int(sections[6])
                hit_end = int(sections[7])
            
                # We can check whether a CDS is completely contained within
                # the region we found.
                for cds in self.transposons[sections[1]].cds_regions:
                    if hit_start <= cds[0] and hit_end >= cds[1]:

                        # We have a hit.  We now need to extract the genomic
                        # positions for the CDS region.
                        # 
                        # The genomic positions are given relative to the strand
                        # on which it occurred.  For negative strand hits then 
                        # the start is higher than the end.
                        #
                        # We are covering an entire CDS
                        genome_start = int(sections[9])
                        genome_end = int(sections[10])
                        chromosome = sections[0]

                        ## TESTING ONLY
                        if chromosome != "chr1":
                            return()


                        strand = sections[8]

                        # Find the region which covers the CDS
                        start_offset = cds[0] - hit_start
                        end_offset = hit_end - cds[1]

                        genome_cds_start = genome_start
                        genome_cds_end = genome_end

                        if strand == "+":
                            genome_cds_start = genome_start + start_offset
                            genome_cds_end = genome_end - end_offset
                        else:
                            genome_cds_start = genome_start - end_offset
                            genome_cds_end = genome_end + start_offset

                        coding_transposon = CodingTransposon(id=sections[1], chromosome=chromosome, cds_length=cds[1]-(cds[0]-1) ,strand=strand, genome_start=genome_start, genome_end=genome_end, genome_cds_start=genome_cds_start, genome_cds_end=genome_cds_end)

                        if not chromosome in self.coding_transposons:
                            self.coding_transposons[chromosome] = []

                        self.coding_transposons[chromosome].append(coding_transposon)

    def extract_cds_sequences(self,genome):
        # Here we go through the putative hits we found before and see if we can
        # find a suitable open reading frame to match the expected location of the
        # cds

        cds_sequences = []

        for chromsome in self.coding_transposons.keys():
            chr_sequence = genome[chromsome]
            for hit in self.coding_transposons[chromsome]:
                # We want to extract the sequence from the 
                # start of the cds to the end of the hit

                potential_orf = chr_sequence[hit.genome_cds_start-1:hit.genome_cds_end]

                # We reverse complement if the hit is on the reverse strand
                if hit.strand == "-":
                    potential_orf = chr_sequence[hit.genome_cds_end-1:hit.genome_cds_start]
                    potential_orf = potential_orf.translate(revcomp_table)[::-1]

                # Now we find ORFs
                orfs = []
                for start_position in [0,1,2]:
                    current_cds_start = start_position
                    current_cds_end = start_position
                    current_position = start_position

                    while (current_position<len(potential_orf)-4):
                        # In the ORF search we're using zero based index positions
                        codon = potential_orf[current_position:current_position+3]
                        if codon in stop_codons:
                            # It's the end of a CDS
                            current_cds_end = current_position-1
                            # We need at least 80% of the full length
                            if ((current_cds_end-(current_cds_start-1))/hit.cds_length) > 0.8:
                                # It's big enough
                                # We use 1 based positions to report
                                orfs.append((current_cds_start+1,current_cds_end+1))
                            
                            current_cds_start = current_position+3
                        
                        current_position += 3

                # See if we have any ORFs
                if not orfs:
                    continue

                # Find the biggest orf
                orfs.sort(key= lambda x: x[1]-x[0])
                biggest_orf = orfs[-1]
                # Convert back to genomic coordinates
                if hit.strand == "+":
                    orf_genomic_start = biggest_orf[0] + hit.genome_cds_start-1
                    orf_genomic_end = biggest_orf[1] + hit.genome_cds_start-1
                else:
                    orf_genomic_start = hit.genome_cds_start - biggest_orf[0]
                    orf_genomic_end = hit.genome_cds_start - biggest_orf[1]

                if not potential_orf[biggest_orf[0]-1:biggest_orf[1]]:
                    breakpoint()

                cds_sequences.append(ExtractedCDS(
                    id=hit.id, 
                    type=self.transposons[hit.id].type,
                    subtype=self.transposons[hit.id].subtype,
                    chromosome=chromsome,
                    strand=hit.strand,
                    repeat_cds_genome_start=hit.genome_cds_start,
                    repeat_cds_genome_end=hit.genome_cds_end,
                    genome_orf_start=orf_genomic_start,
                    genome_orf_end=orf_genomic_end,
                    cds_sequence=potential_orf[biggest_orf[0]-1:biggest_orf[1]]
                ))
        return cds_sequences

class ExtractedCDS:
    def __init__ (self, id, type, subtype, chromosome, strand, repeat_cds_genome_start, repeat_cds_genome_end, genome_orf_start, genome_orf_end, cds_sequence):
        self.id = id
        self.type = type
        self.subtype =subtype
        self.chromosome = chromosome
        self.strand = strand
        self.repeat_cds_genome_start = repeat_cds_genome_start
        self.repeat_cds_genome_end = repeat_cds_genome_end
        self.genome_orf_start = genome_orf_start
        self.genome_orf_end = genome_orf_end
        self.cds_sequence = cds_sequence

class CodingTransposon:

    def __init__ (self, id, chromosome, strand, cds_length, genome_start, genome_end, genome_cds_start, genome_cds_end):
        self.id = id
        self.chromosome = chromosome
        self.strand = strand
        self.cds_length = cds_length
        self.genome_start = genome_start
        self.genome_end = genome_end
        self.genome_cds_start = genome_cds_start
        self.genome_cds_end = genome_cds_end