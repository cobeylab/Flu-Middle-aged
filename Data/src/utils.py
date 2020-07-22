import calc_similarities

def get_seg_sites(viruses):
	WT = ''
	for virus in viruses:
		if virus.split("\n")[0].find("WT") >= 0: 
			WT = virus.split("\n")[1]
			break
			
	seg_sites = set()
	for v in range(len(viruses)):
		seq = viruses[v].split("\n")[1]
		for s in range(len(WT)):
			if WT[s] != seq[s]:
				seg_sites.add(s+1)
	seg_sites = sorted(list(seg_sites))
	seg_sites = [[i] for i in seg_sites]
	
	wtAA = []
	for s in seg_sites:
		wtAA.append(WT[s[0]-1])
	
	return seg_sites, wtAA, WT

def make_name_for_imprinting_to_aa(seg_sites, pngsLocs, epitopes):
	names = []

	for i in range(len(epitopes)):
		seg_sites_1 = [str(x[0]) for x in seg_sites if x[0] in epitopes[i]]
		pngs_1 =      [x[0]+str(x[1]) for x in pngsLocs if x[1] in epitopes[i]]
		imprinting_1 = [ 'impMut'+x+"_"+y for x in seg_sites_1 for y in pngs_1 ]
	
		names += imprinting_1
		
	names = ','.join(names)
	return (names)
	
def make_name_for_imprinting_to_site_by_pngs(pngsLocs, epitopes):
	names = []

	for i in range(len(epitopes)):
		pngs_1 =      [x[0]+str(x[1]) for x in pngsLocs if x[1] in epitopes[i]]
		imprinting_1 = [ 'impSite_'+x for x in pngs_1 ]
			
		names += imprinting_1
		
	names = ','.join(names)
	return (names)	
	
def make_name_for_imprinting_to_site_by_pngs_seg(seg_sites, pngsLocs, epitopes):
	names = []

	for i in range(len(epitopes)):
		seg_sites_1 = [str(x[0]) for x in seg_sites if x[0] in epitopes[i]]
		pngs_1 =      [x[0]+str(x[1]) for x in pngsLocs if x[1] in epitopes[i]]
		imprinting_1 = [ 'impSite_seg_'+x for x in pngs_1 ]
	
		if len(seg_sites_1) == 0:
			continue
			
		names += imprinting_1
		
	names = ','.join(names)
	return (names)	
	
def make_name_for_imprinting_to_aa_non(seg_sites, epitopes):
	epitopes2 = [j for i in epitopes for j in i] #unlist nested list
	seg_sites_non = [x for x in seg_sites if x[0] not in epitopes2]
	names = ['impMut'+str(x[0]) for x in seg_sites if x[0] not in epitopes2]
	names = ','.join(names)
	
	return (seg_sites_non, names)
	
def calc_many_similarities(y, virus, seqByYear, sites_for_similarity, p, max_age, beginY, endY):
	
	similarities = []
	for site in sites_for_similarity:
		similarity = calc_similarities.calc_Svy(y, virus, seqByYear, site, p, max_age, beginY, endY)
		similarities.append(similarity)
		
	return similarities
	
def calc_many_similarities_conserved(y, virus, seqByYear, sites_for_similarity, p, max_age, beginY, endY, seg_sites):
	
	similarities = []
	seg_sites = [j for i in seg_sites for j in i]
	
	for site in sites_for_similarity:
		site = [x for x in site if x not in seg_sites]
		similarity = calc_similarities.calc_Svy(y, virus, seqByYear, site, p, max_age, beginY, endY)
		similarities.append(similarity)
		
	return similarities		
	
def calc_monosite_similarities(y, virus, seqByYear, seg_sites, p, max_age, beginY, endY, epitopes):
	
	seg_sites = [j for i in seg_sites for j in i]#unlist nested list 
	monosites = [j for i in epitopes for j in i] #unlist nested list 
	monosites = [x for x in monosites if x not in seg_sites]
	
	similarity = calc_similarities.calc_Svy(y, virus, seqByYear, monosites, p, max_age, beginY, endY)
	
	return similarity
	
def calc_epitope_similarities(y, virus, seqByYear, seg_sites, p, max_age, beginY, endY, epitopes):
	
	unlisted_epitope = [j for i in epitopes for j in i] #unlist nested list 
	similarity = calc_similarities.calc_Svy(y, virus, seqByYear, unlisted_epitope, p, max_age, beginY, endY)
	
	return similarity
	
	
def calc_many_imprintings(y, wtAA, seqByYear, seg_sites, p, max_age, beginY, endY):
	
	imprintings = []
	for site in seg_sites:
		imprinting = calc_similarities.calc_Ivy(y, wtAA, seqByYear, site, p, max_age, beginY, endY)
		imprintings.append(imprinting)
		
	return imprintings
	
def add_to_line(vars):
	oneline = ''
	for i in range(len(vars)):
		oneline += str(vars[i]) + ","
	return oneline

def call_many_simGlys(y, virus, seqByYear, aaLocs, pngsLocs, p, max_age, delta_opt, beginY, endY, epitopes):
	simGlys = []

	for epitope in epitopes:
		aaLoc = [x[0] for x in aaLocs if x[0] in epitope]
		pngsLoc = [x[1] for x in pngsLocs if x[1] in epitope]

		for aa1 in aaLoc:
			for pngs1 in pngsLoc:

				simGly = calc_similarities.calc_dvySvy(y, virus, seqByYear, [aa1], [pngs1], p, max_age, delta_opt, beginY, endY)
				simGlys.append(simGly)
			
	return (simGlys)
	
def call_many_simGlys_aggregated(y, virus, seqByYear, aaLocs, pngsLocs, p, max_age, delta_opt, beginY, endY, epitopes):
	simGlys = []

	for epitope in epitopes:
		aaLoc = [x[0] for x in aaLocs if x[0] in epitope]
		pngsLoc = [x[1] for x in pngsLocs if x[1] in epitope]

		if len(aaLoc) == 0:
			continue
			
		for pngs1 in pngsLoc:
			
			simGly = calc_similarities.calc_dvySvy(y, virus, seqByYear, aaLoc, [pngs1], p, max_age, delta_opt, beginY, endY)			
			simGlys.append(simGly)
			
	return (simGlys)	
	
def call_many_similarities_aggregated(y, virus, seqByYear, aaLocs, p, max_age, beginY, endY, epitopes):
	simGlys = []

	for epitope in epitopes:
		aaLoc = [x[0] for x in aaLocs if x[0] in epitope]

		if len(aaLoc) == 0:
			continue
			                       
		simGly = calc_similarities.calc_Svy(y, virus, seqByYear, aaLoc, p, max_age, beginY, endY)			
		simGlys.append(simGly)
			
	return (simGlys)	

	
def call_many_simGlys_by_pngs(y, virus, seqByYear, pngsLocs, p, max_age, delta_opt, beginY, endY, epitopes):
	simGlys = []
	
	for epitope in epitopes:

		pngsLoc = [x[1] for x in pngsLocs if x[1] in epitope]

		for pngs1 in pngsLoc:
			
			simGly = calc_similarities.calc_dvySvy(y, virus, seqByYear, epitope, [pngs1], p, max_age, delta_opt, beginY, endY)			
			simGlys.append(simGly)
			
	return (simGlys)
	
def call_many_simGlys_by_pngs_conserved(y, virus, seqByYear, pngsLocs, p, max_age, delta_opt, beginY, endY, epitopes, seg_sites):
	simGlys = []
	seg_sites = [j for i in seg_sites for j in i]#unlist nested list 
	
	for epitope in epitopes:

		pngsLoc = [x[1] for x in pngsLocs if x[1] in epitope]
		epitope = [x for x in epitope if x not in seg_sites]

		for pngs1 in pngsLoc:
			
			simGly = calc_similarities.calc_dvySvy(y, virus, seqByYear, epitope, [pngs1], p, max_age, delta_opt, beginY, endY)			
			simGlys.append(simGly)
			
	return (simGlys)

def call_many_glyProbs(y, virus, seqByYear, pngsLocs, p, max_age, delta_opt, beginY, endY):
	
	glyProbs = []
	for pngsLoc in pngsLocs:
		glyProb = calc_similarities.calc_GEvy(y, virus, seqByYear, [pngsLoc[1]], p, max_age, delta_opt, beginY, endY)
		glyProbs.append(glyProb)
	
	return glyProbs
	
	
	
	
	
	