import csv
import os
import re
import subprocess
import pysam
import glob
################################################
# You can use this code and put it in your own script
class ParseFastQ(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for rec in parser:
            ... do something with rec ...
 
        rec is tuple: (seqHeader,seqStr,qualHeader,qualStr)
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'r')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
         
    def __iter__(self):
        return self
     
    def __next__(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: tuple: (seqHeader,seqStr,qualHeader,qualStr)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else: 
                elemList.append(None)
         
        # ++++ Check Lines For Expected Form ++++
        trues = [bool(x) for x in elemList].count(True)
        nones = elemList.count(None)
        # -- Check for acceptable end of file --
        if nones == 4:
            raise StopIteration
        # -- Make sure we got 4 full lines of data --
        assert trues == 4,\
               "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        # -- Make sure we are in the correct "register" --
        assert elemList[0].startswith(self._hdSyms[0]),\
               "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        assert elemList[2].startswith(self._hdSyms[1]),\
               "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        # -- Make sure the seq line and qual line have equal lengths --
        assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
               Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
         
        # ++++ Return fatsQ data as tuple ++++
        return tuple(elemList)
##########################################################################
def pileup(path):
    #test file, replaced with the sorted.bam you are using. Make sure it is indexed! (Use samtools index yourbam.sorted.bam)
    samfile = pysam.AlignmentFile(path, "rb")

    #Since our reference only has a single sequence, we're going to pile up ALL of the reads. Usually you would do it in a specific region (such as chromosome 1, position 1023 to 1050 for example)
    for pileupcolumn in samfile.pileup():
        #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        #use a dictionary to count up the bases at each position
        ntdict = {}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                # You can uncomment the below line to see what is happening in the pileup. 
                #print ('\tbase in read %s = %s' % (pileupread.alignment.query_name, pileupread.alignment.query_sequence[pileupread.query_position]))
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                ########## ADD ADDITIONAL CODE HERE ############# 
                # Populate the ntdict with the counts of each base 
                # This dictionary will hold all of the base read counts per nucletoide per position.
                try:
                    ntdict[pileupread.alignment.query_sequence[pileupread.query_position]] += 1
                except KeyError:
                    ntdict[pileupread.alignment.query_sequence[pileupread.query_position]] = 1
                # Use the dictionary to calculate the frequency of each site, and report it if if the frequency is NOT  100% / 0%. 
                #############################################
        #Creates returnList that has: number of reads, percentage of mutation, position of mutation, mutation type
        baseCount = pileupcolumn.n;
        for base in ntdict:
            freq = round((ntdict[base]/baseCount)*100)
            if(freq != 0 and freq != 100):
                returnList = [baseCount,round((ntdict[min(ntdict,key=ntdict.get)]/baseCount)*100),pileupcolumn.pos,min(ntdict,key=ntdict.get)]
                return returnList
                #print(base + " frequency is: "+str(freq)+" at position "+str(pileupcolumn.pos))
        #print (ntdict)
    samfile.close()
##########################################################################
#Makes fastqs directory
os.makedirs('fastqs', exist_ok = True)

#Create and populate a dictionary with patient names and barcode
patientDict = {}
with open('harrington_clinical_data.txt', 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    next(file)
    for row in reader:
        patientDict.update({row[0]:row[2]})

#trimBeg removes the barcode and associated quality scores from the fastqObj
#Input: barcode and fastqObj
#Output: fastqObj
def trimBeg(barcode, fastqObj):
    fastqList = [fastqObj[0],fastqObj[1].removeprefix(barcode),fastqObj[2],fastqObj[3][len(barcode):]]
    return fastqList

#trimEnd trims off the end of the quality scores given consecutive D or F quality score, and associated sequence bases
#Input: fastqObj 
#Output: fastqObj
def trimEnd(fastqObj):
    fastqList = [fastqObj[0],fastqObj[1],fastqObj[2],fastqObj[3]]
    fastqList[3] = re.split("[D,F][D,F]",fastqList[3])[0]
    fastqList[1] = fastqList[1][0:len(fastqList[3])]
    return fastqList

#Generates trimmed fastq files for each patient given the pooled sequences fastq
for patient in patientDict:
    fastqFile = ParseFastQ("./hawkins_pooled_sequences.fastq")
    barcode = patientDict[patient]
    for fastqObj in fastqFile:
        if(fastqObj[1].startswith(barcode)):
            fastqList = trimEnd(trimBeg(barcode, fastqObj))
            with open("./fastqs/"+patient+"_trimmed.fastq", 'a') as file:
                for line in fastqList:
                    file.write(line + '\n')

#Creates directory for SAM files, creates list of associated fastq files, and indexs references FASTA for BWA
os.makedirs('./sams', exist_ok = True)
patientList = os.listdir('./fastqs')
subprocess.run(["bwa","index","dgorgon_reference.fa"])

#Creates SAM files
for patientFile in patientList:
    cmd = "bwa mem dgorgon_reference.fa ./fastqs/"+patientFile+" > ./sams/"+patientFile.partition("_")[0]+".sam"
    subprocess.run(cmd,shell=True)

#Creates directory for unsorted BAM files, and creates list of associated SAM files
os.makedirs('./bam', exist_ok = True)
samList = os.listdir('./sams')

#Creates unsorted BAM files
for samFile in samList:
    cmd = "samtools view -bS ./sams/"+samFile+" > ./bam/"+samFile.partition(".")[0]+".bam"
    subprocess.run(cmd,shell=True)

#Creates directory for sorted BAM files, and creates list of associated unsorted BAM files
os.makedirs('./bams', exist_ok = True)
bamList = os.listdir('./bam')

#Creates sorted BAM files
for bamFile in bamList:
    cmd = "samtools sort -m 100M -o ./bams/"+bamFile.partition(".")[0]+".sorted.bam"+" ./bam/"+bamFile
    subprocess.run(cmd,shell=True)
    
    cmd = "samtools index ./bams/"+bamFile.partition(".")[0]+".sorted.bam"
    subprocess.run(cmd,shell=True)

#Removes directories with SAM files and unsorted BAM files
cmd = "rm -r ./sams"
subprocess.run(cmd,shell=True)
cmd = "rm -r ./bam"
subprocess.run(cmd,shell=True)

#Create and populate a dictionary with patient names and mold color
patientDict = {}
with open('harrington_clinical_data.txt', 'r') as file:
    reader = csv.reader(file, delimiter='\t')
    next(file)
    for row in reader:
        patientDict.update({row[0]:row[1]})

#Generates report with SNP info for each patient
bamList = glob.glob('./bams/*.sorted.bam')
report = open('report.txt', 'a')
for bamFile in bamList:
    variantList = pileup(bamFile)
    name = os.path.basename(bamFile).partition(".")[0]
    report.write("Sample "+str(name)+" had a "+str(patientDict[name])+" mold, "+str(variantList[0])+" reads, and "+str(variantList[1])+"% of the reads at position "+str(variantList[2])+" had the mutation "+str(variantList[3])+".\n")
report.close()
