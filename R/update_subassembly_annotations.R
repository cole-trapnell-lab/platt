udpate_annotations = function(cds){
  colData(cds)$updated_cell_type = case_when(
    colData(cds)$subassembly_group %in% c("39-5") ~ "delta cell", # sst2
    colData(cds)$subassembly_group %in% c("39-2") ~ "alpha cell", #gcg
    colData(cds)$subassembly_group %in% c("39-1") ~ "beta cell", #ins
    #colData(notochord_cds)$subassembly_group %in% c("39-6") ~ "ventral habenular nucleus", # just one paper (schier) says this, but they looked at 7dpf
    colData(cds)$subassembly_group %in% c("39-4") ~ "hypophysis", #prl
    colData(cds)$subassembly_group %in% c("39-3", "39-8", "39-9") ~ "adenohypophysis", # pomca

    # Ear: mostly from https://elifesciences.org/articles/82978#fig3
    ## Ear support cells:
    colData(cds)$subassembly_group %in% c("36-6", "36-1") ~ "early otic vesicle", # pax5, cldna
    colData(cds)$subassembly_group %in% c("36-7") ~ "otic nonsensory epithelium", # dachc
    colData(cds)$subassembly_group %in% c("36-5") ~ "saccular macula support cell", # tectb/tecta, otog, missing pax5
    colData(cds)$subassembly_group %in% c("36-8") ~ "utricular macula support cell", # tectb/tecta, otog + pax5 (see https://anatomypubs.onlinelibrary.wiley.com/doi/full/10.1002/dvdy.10073)
    colData(cds)$subassembly_group %in% c("36-10") ~ "crista support cell", # tectb, zlpd1a
    colData(cds)$subassembly_group %in% c("36-4") ~ "semicircular canal", # chsy1 https://www.cell.com/ajhg/fulltext/S0002-9297(10)00519-7
    colData(cds)$subassembly_group %in% c("36-11") ~ "Unknown otic support cell",
    colData(cds)$subassembly_group %in% c("36-3", "36-9", "36-2") ~ "undifferentiated otic support cell",
    ## Ear hair cells:
    colData(cds)$subassembly_group %in% c("44-3") ~ "extrastriolar hair cell", #rtn4r, ush1c
    colData(cds)$subassembly_group %in% c("44-4") ~ "hair cell/support cell doublets?", #mix of many shi genes, small cluster
    colData(cds)$subassembly_group %in% c("44-1") ~ "striolar hair cell", # pvalb9, cabp2b
    colData(cds)$subassembly_group %in% c("44-2") ~ "crista hair cell", # rbm24a (https://academic.oup.com/cardiovascres/article/94/3/418/354969), pchd15a (https://zfin.org/ZDB-PUB-140813-4)

    ## Roof plate, subcomissural organ:
    colData(cds)$subassembly_group %in% c("42-8") ~ "subcommissural organ", #scospondin"
    #colData(cds)$subassembly_group %in% c("42-7", "42-1", "42-1") ~ "Midbrain roof plate" # gpc5a
    #colData(cds)$subassembly_group %in% c() ~ "rhombomere" # otx2b

    ## Cranial & pec fin muscle
    #colData(cds)$subassembly_group %in% c("17-4") ~ "mandibular muscle", #dlx4b -

    colData(cds)$subassembly_group %in% c("17-10") ~ "pectoral fin muscle", #lmx1a, myhz2, myhz1.1, mylz3  (see https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001168,  https://www.molbiolcell.org/doi/10.1091/mbc.e13-08-0486 and Rauch et al., 2003

    ## Neural crest (based on L&S markers + Soldatov et al 2019):
    colData(cds)$subassembly_group %in% c("6-10", "6-1") ~ "neural crest (pre-migratory)", # dlx5a, arhgef7a
    colData(cds)$subassembly_group %in% c("18-2") ~ "neural crest (migrating, multi-potent)", #
    colData(cds)$subassembly_group %in% c("18-4", "18-5", "13-6", "18-10") ~ "neural crest (pigment fated)",
    colData(cds)$subassembly_group %in% c("13-9") ~ "xanthophore", #bscl2l

    # see https://elifesciences.org/articles/60005 for lots of good markers of sox10+ cells
    colData(cds)$subassembly_group %in% c("6-26", "6-9", "6-35") ~ "enteric progenitors", # ngrfb, ret, gfra1a, gfra3
    #colData(cds)$subassembly_group %in% c() ~ "autonomic neuron progenitors (neural crest-derived)", # hand2, lmo1 plus lamc3, tnn (axon guidance)


    colData(cds)$subassembly_group %in% c("18-7", "18-8", "18-3", "18-1") ~ "neural crest (sensory neuron fated)", # ngrfb
    colData(cds)$subassembly_group %in% c("6-18", "6-3", "6-2", "6-22", "6-32", "6-8", "6-33") ~ "neural crest (mesenchyme fated)", # meox1, foxc1a (used to be called "pharyngeal arch (contains muscle, early cartilage)"
    colData(cds)$subassembly_group %in% c("6-39", "6-12", "6-36") ~ "chondrocytes (pharyngeal-arch derived)", # meox1, foxc1a
    colData(cds)$subassembly_group %in% c("6-42", "6-25", "6-6") ~ "osteocytes (pharyngeal-arch derived)", # panx3
    colData(cds)$subassembly_group %in% c("6-19") ~ "endothelium (pharyngeal-arch derived)", #
    colData(cds)$subassembly_group %in% c("6-33") ~ "tenocyte (pharyngeal-arch derived)", # tnmd, cdon

    #colData(cds)$subassembly_group %in% c("") ~ "doublets", # probably PA + endothelium

    TRUE ~ colData(cds)$cell_type_sub
  )
  return(cds)
}
