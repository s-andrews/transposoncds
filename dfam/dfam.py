#!/usr/bin/env python
import gzip

class Dfam:

    def __init__(self,file):
        self.file = file


    def parseFile(self):
        # We need to parse transposons which have CDS features in them.  All
        # others are irrelevant.  The only information we need from them is
        # 
        # From: ID   DF0000018; SV 4; linear; DNA; STD; UNC; 2781 BP.
        #
        # 1. Their Dfam ID and version DF0000018.4
        # 2. Their length
        #
        # From: NM   Charlie1
        #
        # 3. Their name
        #
        # From: KW   DNA/hAT-Charlie.
        # 
        # 4. Their type (DNA)
        # 5. Their subtype (hAT-Charlie)
        #
        # From: FT   CDS             597..2538
        #
        # 6. The start position (597)
        # 7. The end position (2538)

        self.transposons = {}

        with gzip.open(self.file) as infh:
            id = None
            length = None
            name = None
            type = None
            subtype = None
            cds_regions = []

            for line in infh:

                # Check for an ID line
                if line.startswith("ID"):
                    # See if we have a complete entry
                    if id is not None and cds_regions:
                        self.transposons[id] = Transposon(id,length,name,type,subtype,cds_regions)

                    # Wipe the held data
                    id = None
                    length = None
                    name = None
                    type = None
                    subtype = None
                    cds_regions = []

                    sections = line[5:].split("; ")
                    id = sections[0]+sections[1][3:]
                    length = int(sections[-1].split()[0])

                elif line.startswith("NM"):
                    name = line[2:].strip()

                elif line.startswith("KW"):
                    sections = line[2:].strip().split("/")
                    type = sections[0]
                    subtype =sections[1]

                elif line.startswith("FT   CDS"):
                    sections = line[8:].strip().split("..")
                    cds_start = int(sections[0])
                    cds_end = int(sections[1])
                    cds_regions.append((cds_start,cds_end))


class Transposon:
    def __init__(self,id,length,name,type,subtype,cds_regions):
        self.id = id
        self.length = length
        self.name = name
        self.type = type
        self.subtype = subtype
        self.cds_regions = cds_regions
