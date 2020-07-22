import random

def read_epitope_shih(epifName):
	sites = ['a','b','c','d','e']
	epitopes = [[] for i in range(len(sites)) ]

	epif = open(epifName, "r")
	for line in epif:
		each = line.split("\n")[0].split(",")
		idx = sites.index(each[1])
		epitopes[idx].append(int(each[0]))

	epif.close()

	return epitopes

def read_epitope_koel(epifName):

	epitope = []

	epif = open(epifName, "rU")
	for line in epif:
		epitope.append(int(line.split("\n")[0]))

	epif.close()

	return epitope
	
def read_subtype(subf_name, beginY, endY, col):
	subf = open(subf_name, "rU")
	p_subtype_dict = dict()
	for line in subf:
		if line.find("YOB") >= 0:
			continue
		else:
			each = line.split(",")
			yob = each[0]
			p_subtype_dict[int(yob)] = each[col]

	return p_subtype_dict

def read_virus(testv_name):
	testvf_name = "../data/"+testv_name+"_AA_std.fas"
	testvf = open(testvf_name, "rU")
	viruses = []
	vyears = []
	for line in testvf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			testv = each[1] + "\n"
			year = each[2].split("/")[0]
		else:
			testv += line.split("\n")[0]
			viruses.append(testv)
			vyears.append(year)
	return viruses, vyears
	
def read_NA_virus(testvf_name):
	testvf = open(testvf_name, "rU")
	viruses = []
	vyears = []
	for line in testvf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			testv = each[1] + "\n"
			year = each[2].split("/")[0]
		else:
			testv += line.split("\n")[0]
			viruses.append(testv)
			vyears.append(year)
	return viruses, vyears

def sample_nSamp(byYear, nSamp):
	samples = []
	for year in byYear:
		if len(year) > nSamp:
			sample = random.sample(year, nSamp)
		else:
			sample = year
		samples.append(sample)

	return samples

def sequences_byYear(firstY, lastY, seqfName):
	byYear = [[] for i in range(firstY, lastY+1)]

	seqf = open(seqfName, "rU")
	for line in seqf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			y = int(each[2].split("/")[0])
		else:
			if y < firstY or y > lastY:
				continue
				
			seq = line.split("\n")[0] 			
			if (seq.find("-") >= 0 and seq.find("-") < 457):
				continue
				
			byYear[y-firstY].append(line.split("\n")[0])
			
	#for y in byYear:
	#	print (len(y))
	return byYear
	
def sequences_bySeason(firstS, lastS, seqfName):
	bySeason = [[] for i in range(firstS, lastS+1)]
	seq_t = 0
	seq_tm1 = 0
	seqf = open(seqfName, "rU")
	for line in seqf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			y = int(each[2].split("/")[0])
			s = ''
			try:
				m = int(each[2].split("/")[1])
				if (m <= 9):
					s = y
					seq_t += 1
				elif (m > 9):
					s = y+1
					seq_tm1 +=1
			except ValueError:
				continue
				
		else:

			if s == '':
				continue
			if s < firstS or s > lastS:
				continue
				
			seq = line.split("\n")[0] 			
			if (seq.find("-") >= 0 and seq.find("-") < 457):
				continue
								
			bySeason[s-firstS].append(line.split("\n")[0])
			
	#print (seq_t, seq_tm1)
	return bySeason

def sequences_year_to_season(seqByYear):
	byYear = [[] for i in range(firstY, lastY+1)]

	seqf = open(seqfName, "rU")
	for line in seqf:
		if line.find(">") >= 0:
			each = line.split("\n")[0].split("|")
			y = int(each[2].split("/")[0])
		else:
			if y < firstY or y > lastY:
				continue
			byYear[y-firstY].append(line.split("\n")[0])
			
	#for y in range(len(seqByYear)):
		
	
	return byYear

def vacSequences_bySeason(vseqf_name, egg_mutations):
	vseqf = open(vseqf_name, "r")
	vacSeqBySeason = []
	
	for line in vseqf:
		if line.find(">") >= 0:
			continue
		else:
			vacSeqBySeason.append(line.split("\n")[0])
	
	print (len(vacSeqBySeason), len(egg_mutations))
	return vacSeqBySeason
	
def read_scales(scalef):
	scale_data = open(scalef, "U")
	
	scales = [[] for i in range(1918, 2016+1)]
	
	for line in scale_data:
		if line.find("birth") >= 0:
			continue
		byear = float(line.split(",")[0])
		season = float(line.split(",")[1])
		scale = float(line.split(",")[2])
		
		if season == 2009.5:
			continue
		
		if byear < 2017:
			scales[int(byear)-1918].append(scale)
		
	return (scales)
	
def read_vac_cov_by_byear(fname):
	vac_cov_f = open(fname, "U")
	
	vac_cov_by_byear = []
	for line in vac_cov_f:
		each = line.split("\n")[0].split(",")
		vac_cov_y = []
		for e in each:
			vac_cov_y.append(float(e))
		vac_cov_by_byear.append(vac_cov_y[1:])
		
	return (vac_cov_by_byear)
		
	
def read_naive_for_vac(is_full):
	
	
	
	if is_full==0:
		naive_for_vac = [0.0 for i in range(0, 20+1) ] 
		naive_for_vac[1] = 1.0
		naive_for_vac[2] = 0.37
		naive_for_vac[3] = 0.07
		naive_for_vac[4] = 0.04
		naive_for_vac[5] = 0.02
		
	elif is_full == 1:
		naive_for_vac = [1.0 for i in range(0, 20+1) ] 

	
	return (naive_for_vac)
	
	
	
def read_na_fasta(naf_name):
	naf = open(naf_name, "U")
	for line in naf:
		if line.find(">") >= 0:
			continue
		else:
			return line.split("\n")[0]
	
def read_vac_cov_for1617 (vacf_name):
	vacf = open(vacf_name, "U")
	vc_dict = dict()
	for line in vacf:
		if line.find("Age") >= 0:
			continue
		each = line.split("\n")[0].split(",")
		age = int(each[0])
		vac_cov = float(each[2])
		vc_dict[age] = vac_cov
	
	return vc_dict
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	