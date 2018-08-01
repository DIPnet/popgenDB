
USETEXTLINKS = 1
STARTALLOPEN = 0
WRAPTEXT = 1
PRESERVESTATE = 0
HIGHLIGHT = 1
ICONPATH = 'file:////Users/eric/github/popgenDB/sims_for_structure_paper/2PopDNAnorec_0.5_1000/'    //change if the gif's folder is a subfolder, for example: 'images/'

foldersTree = gFld("<i>ARLEQUIN RESULTS (2PopDNAnorec_1_76.arp)</i>", "")
insDoc(foldersTree, gLnk("R", "Arlequin log file", "Arlequin_log.txt"))
	aux1 = insFld(foldersTree, gFld("Run of 31/07/18 at 17:02:31", "2PopDNAnorec_1_76.xml#31_07_18at17_02_31"))
	insDoc(aux1, gLnk("R", "Settings", "2PopDNAnorec_1_76.xml#31_07_18at17_02_31_run_information"))
		aux2 = insFld(aux1, gFld("Genetic structure (samp=pop)", "2PopDNAnorec_1_76.xml#31_07_18at17_02_31_pop_gen_struct"))
		insDoc(aux2, gLnk("R", "AMOVA", "2PopDNAnorec_1_76.xml#31_07_18at17_02_31_pop_amova"))
		insDoc(aux2, gLnk("R", "Pairwise distances", "2PopDNAnorec_1_76.xml#31_07_18at17_02_31_pop_pairw_diff"))
