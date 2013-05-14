/*=============================================================================

    SPIMAP - Species Informed Max A Posteriori Phylogenetic Reconstruction

    Matt Rasmussen
    Copyright 2010-2011

    modified for WGD by Charles-Elie Rabier
    2011-2013

=============================================================================*/

// c++ headers
#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>
#include <memory>

// third party headers
#include <gsl/gsl_errno.h>

// spidir headers
#include "common.h"
#include "ConfigParam.h"
#include "logging.h"
#include "model.h"
#include "model_params.h"
#include "newick.h"
#include "parsimony.h"
#include "parsing.h"
#include "phylogeny.h"
#include "search.h"
#include "seq.h"
#include "seq_likelihood.h"
#include "Sequences.h"
#include "treevis.h"
#include "WGD.h"

#define VERSION_TEXT "1.2"
#define VERSION_INFO  "\
SPIMAP " VERSION_TEXT " \n\
SPecies Informed Max A Posteriori gene tree reconstruction \n\
Matt Rasmussen 2010-2011\n\
CSAIL, MIT \n\
\n\
Modified by Charles-Elie Rabier (2011) for Whole Genome Duplications\n\
Citation:\n\
A Bayesian Approach for Fast and Accurate Gene Tree Reconstruction.\n\
Matthew D. Rasmussen and Manolis Kellis.\n\
Molecular Biology and Evolution. 2010.\n\
"


//=============================================================================

using namespace std;
using namespace spidir;

// debug options level
const int DEBUG_OPT = 1;


// parsing command-line options
class SpidirConfig
{
public:

    SpidirConfig() 
    {
        makeParser();
    }

    void makeParser()
    {
        config.clear();

        // input/output
	config.add(new ConfigParam<string>
		   ("-a", "--align", "<alignment fasta>", &alignfile, 
		    "sequence alignment in fasta format"));
	config.add(new ConfigParam<string>
		   ("-S", "--smap", "<species map>", &smapfile, 
		    "gene to species map"));
	config.add(new ConfigParam<string>
		   ("-s", "--stree", "<species tree>", &streefile, 
		    "species tree file in newick format"));
	config.add(new ConfigParam<string>
		   ("-p", "--param", "<params file>", &paramsfile, 
		    "substitution rate parameters file"));
	config.add(new ConfigParam<string>
		   ("-o", "--output", "<output filename prefix>", 
		    &outprefix, "spimap",
		    "prefix for all output filenames"));
	config.add(new ConfigSwitch
		   ("-r", "--recon", 
		    &outputRecon,
		    "Output reconciliation"));
    

        // sequence model
	config.add(new ConfigParamComment("Sequence evolution model"));
	config.add(new ConfigParam<float>
		   ("-k", "--kappa", "<transition/transversion ratio>", 
		    &kappa, -1.0,
		    "used for HKY model (default=estimate)"));
	config.add(new ConfigParam<string>
		   ("-f", "--bgfreq", "<A freq>,<C ferq>,<G freq>,<T freq>", 
		    &bgfreqstr, "",
		    "background frequencies (default: estimate)"));

        // dup/loss model
	config.add(new ConfigParamComment("Dup/loss evolution model"));
	config.add(new ConfigParam<float>
		   ("-D", "--duprate", "<duplication rate>", 
		    &duprate, 0.1,
		    "rate of a gene duplication (default=0.1)"));
	config.add(new ConfigParam<float>
		   ("-L", "--lossrate", "<loss rate>", 
		    &lossrate, 0.1,
		    "probability of loss (default=0.1)"));
	config.add(new ConfigParam<float>
		   ("-LR", "--lineagesatroot", "<param for lineages at the root>                     ",&q, 0.5,
		    "parameter for the geometric law, concerning the number of lineages at the root (default=0.5)"));
	config.add(new ConfigParam<int>
		    ("-OB","--observingsomething", "<obssmthg>", &observingsomething, 1,
		    "1 for conditionning on observing something (default: 1) or 2 for at least one on the left and at least one on the right"));
	config.add(new ConfigParam<float>
		   ("-P", "--pretime", "<pre-speciation time parameter>", 
		    &pretime, 1.00,
		    "lambda param of pre-speciation distribution (default=1.0)"));

        // search options
	config.add(new ConfigParamComment("Search"));
	config.add(new ConfigParam<int>
		   ("-i", "--niter", "<# iterations>", 
		    &niter, 100, 
		    "number of iterations"));
	config.add(new ConfigParam<int>
		   ("", "--quickiter", "<quick iterations>", 
		    &quickiter, 50,
		    "number of subproposals (default=50)"));
	config.add(new ConfigParam<int>
		   ("-b", "--boot", "<# bootstraps>", 
		    &bootiter, 1,
		    "number of bootstraps to perform (default: 1)"));
        config.add(new ConfigParam<int>
		   ("-g", "--proposal-gene-topology", "<proposal type for gene tree topology>", 
		    &propid, 1,
		    "1 for spr-neighbor (default: 1), 0 for nni, 2 for SubtreeSlide"));
	 config.add(new ConfigParam<int>
		    ("","--mcmc", "<mcmc>", 
		    &method, 1,
		    "1 for MCMC (default: 1) or 0 for MAP"));
    
        // misc
	config.add(new ConfigParamComment("Miscellaneous", DEBUG_OPT));
	config.add(new ConfigParam<string>
		   ("", "--search", "default|sprnbr|nosprnbr", 
		    &search, "default", 
		    "search algorithm", DEBUG_OPT));
	config.add(new ConfigParam<string>
		   ("", "--prior", "spimap|none", 
		    &prioropt, "spimap",
		    "function for prior (default=spimap)", DEBUG_OPT));
	config.add(new ConfigParam<string>
		   ("-c", "--correct", "<correct tree file>", &correctFile, ""
		    "check if correct tree is visited in search",
                    DEBUG_OPT));
	config.add(new ConfigParam<int>
		   ("", "--prior_samples", "<number of samples>",
		    &priorSamples, 100,
		    "number of samples to use in branch prior integration (default: 100)", DEBUG_OPT));
	config.add(new ConfigSwitch
		   ("", "--prior_exact", 
		    &priorExact,
		    "Use an exact calculation of branch prior", DEBUG_OPT));
	config.add(new ConfigParam<int>
		   ("", "--lkiter", "<max number of likelihood iterations>",
		    &lkiter, 10,
		    "max number of iterations in maximum likelihood fitting (default: 10)", DEBUG_OPT));
        config.add(new ConfigParam<float>
                   ("", "--minlen", "<length>",
                    &minlen, 0.0,
                    "minimum branch length allowed", DEBUG_OPT));
        config.add(new ConfigParam<float>
                   ("", "--maxlen", "<length>",
                    &maxlen, 10.0,
                    "maximum branch length allowed", DEBUG_OPT));
        config.add(new ConfigParam<int>
		   ("-x", "--seed", "<random seed number>", 
		    &seed, 0, 
		    "use 0 to use time as seed", DEBUG_OPT));

        // help information
	config.add(new ConfigParamComment("Information"));
	config.add(new ConfigParam<int>
		   ("-V", "--verbose", "<verbosity level>", 
		    &verbose, LOG_LOW, 
		    "verbosity level 0=quiet, 1=low, 2=medium, 3=high"));
	config.add(new ConfigParam<string>
		   ("", "--log", "<log filename>", &logfile, "", 
		    "log filename.  Use '-' to display on stdout."));
	config.add(new ConfigSwitch
		   ("", "--treeSampled", 
		    &keepTreeSampled,
		    "Output treeSampled"));
	config.add(new ConfigSwitch
		   ("", "--informationduploss", 
		    &keepDupLoss,
		    "Output losses, losses at WGD, duplications, duplications at WGD"));
	config.add(new ConfigSwitch
		   ("-v", "--version", &version, "display version information"));
	config.add(new ConfigSwitch
		   ("-h", "--help", &help, 
		    "display help information"));
	config.add(new ConfigSwitch
		   ("", "--help-debug", &help_debug, 
		    "display help information about debug options"));        
    }

    int parseArgs(int argc, char **argv)
    {
	// parse arguments
	if (!config.parse(argc, (const char**) argv)) {
	    if (argc < 2)
		config.printHelp();
	    return 1;
	}

	// display help
	if (help_debug) {
	    config.printHelp(stderr, DEBUG_OPT);
	    return 1;
	}
    
	// display help
	if (help) {
	    config.printHelp();
	    return 1;
	}
    
	// display version info
	if (version) {
	    printf(VERSION_INFO);
	    return 1;
	}
    	

	if (duprate == lossrate)
	    lossrate *= .98;

	return 0;
    }


  void printarguments()
  {

    printLog(LOG_LOW, "SPIMAP executed with the following parameters\n");
    printLog(LOG_LOW, "-propGT (0 for NNI and 1 for SPR) %d\n", propid);
    printLog(LOG_LOW, "-a %s\n", alignfile.c_str());
    printLog(LOG_LOW, "-S %s\n", smapfile.c_str());
    printLog(LOG_LOW, "-s %s\n", streefile.c_str());
    printLog(LOG_LOW, "-p %s\n", paramsfile.c_str());
    printLog(LOG_LOW, "-o %s\n", outprefix.c_str());
    printLog(LOG_LOW, "-r (1 true, 0 false) %d\n", outputRecon);
    printLog(LOG_LOW, "-k  %f\n", kappa);
    printLog(LOG_LOW, " -f %s\n", bgfreqstr.c_str());
    printLog(LOG_LOW, "-duprate %f\n", duprate);
    printLog(LOG_LOW, "-lossrate %f\n", lossrate);  
    printLog(LOG_LOW, "--lineagesatroot %f\n", q);
    printLog(LOG_LOW, "--observingsomething (1 for observingsomething and 2 for at least one on the left and at least one on the right %d\n", observingsomething);
    printLog(LOG_LOW, "-P %f\n", pretime);
    printLog(LOG_LOW, "-niter %d\n", niter);
    printLog(LOG_LOW, "-- quickiter %d\n", quickiter);
    printLog(LOG_LOW, "-b %d\n", bootiter);
    printLog(LOG_LOW, "-g (0 for NNI, 1 for SPR, 2 for SubtreeSlide) %d\n", propid);
    printLog(LOG_LOW, "--mcmc (1 for MCMC and 0 for MAP) %d\n", method);
    printLog(LOG_LOW, "-x %d\n", seed);
    printLog(LOG_LOW, "comment %s\n", search.c_str());
    printLog(LOG_LOW, "-c %s\n", correctFile.c_str());
    printLog(LOG_LOW, "--prior %s\n", prioropt.c_str());
    printLog(LOG_LOW, "--prior_samples %d\n", priorSamples);
    printLog(LOG_LOW, "--prior_exact %d\n",  priorExact);
    printLog(LOG_LOW, "--lkiter %d\n",  lkiter);
    printLog(LOG_LOW, "--minlen %f\n", minlen);
    printLog(LOG_LOW, "--maxlen %f\n", maxlen);
    printLog(LOG_LOW, "-V %d\n", verbose);
    printLog(LOG_LOW, "--treeSampled (1 true, 0 false) %d\n", keepTreeSampled);
    printLog(LOG_LOW, "--informationduploss (1 true, 0 false) %d\n", keepDupLoss);
    printLog(LOG_LOW, "-v %d\n", version);
    printLog(LOG_LOW, "-h %d\n", help);
    printLog(LOG_LOW, "--help-debug %d\n", help_debug);
    printLog(LOG_LOW, "--log %s\n", logfile.c_str());
    printLog(LOG_LOW, "\n\n");
  }
  


    ConfigParser config;

    // input/output
    string alignfile;
    string smapfile;
    string streefile;
    string paramsfile;
    string outprefix;
    bool outputRecon;

    // sequence model
    float kappa;
    string bgfreqstr;

    // dup/loss model
    float duprate;
    float lossrate;
    float pretime;
    float q;
    int observingsomething;
    
    // search
    int niter;
    int quickiter;
    int bootiter;
    int propid;
    int method;

    // misc
    int seed;
    string search;
    string correctFile;
    string prioropt;
    int priorSamples;
    bool priorExact;
    int lkiter;
    float minlen;
    float maxlen;

    // help/information
    int verbose;
    bool keepTreeSampled;
    bool keepDupLoss;
    bool version;
    bool help;
    bool help_debug;
    string logfile;
 
};


// perform bootstrapping
bool bootstrap(Sequences *aln, string *genes, TreeSearch *search,
	       int bootiter, string outprefix)
{
    // bootstrap
    if (bootiter > 1) {  
	Tree *boottree;
    
	string bootFilename = outprefix + ".boot.trees";
	string bootAlignFilename = outprefix + ".boot.align";
	FILE *bootfile = NULL;
	FILE *bootAlignfile = NULL;
            
	if (! (bootfile = fopen(bootFilename.c_str(), "w"))) {
	    printError("cannot open '%s' for writing", bootFilename.c_str());
	    return false;
	}
            
	if (! (bootAlignfile = fopen(bootAlignFilename.c_str(), "w"))) {
	    printError("cannot open '%s' for writing", bootAlignFilename.c_str());
	    return false;
	}

	// create blank alignment for bootstrapping
	Sequences aln2;
	aln2.alloc(aln->nseqs, aln->seqlen);
	for (int i=0; i<aln->nseqs; i++)
	    aln2.names[i] = aln->names[i];

	for (int i=1; i<=bootiter; i++) {
	    printLog(LOG_LOW, "bootstrap %d of %d\n", i, bootiter);
	    resampleAlign(aln, &aln2);
            
	    boottree = search->search(NULL, genes, 
				      aln2.nseqs, aln2.seqlen, aln2.seqs);

	    boottree->setLeafNames(genes);
	    writeNewickTree(bootfile, boottree, 0, true);
	    fprintf(bootfile, "\n");
	    fflush(bootfile);
	    delete boottree;            

	    // DEBUG
	    writeFasta(bootAlignfile, &aln2);
	}

	fclose(bootfile);
	fclose(bootAlignfile);
    } 
    
    return true;
}



int main(int argc, char **argv)
{
    SpidirConfig c;
    int ret = c.parseArgs(argc, argv);
    if (ret)
	return ret;

    //=======================================================
    // setup gsl
    gsl_set_error_handler_off();

    
    //=======================================================
    // logging
    
    // use default log filename
    if (c.logfile == "")
        c.logfile = c.outprefix + ".log";
    
    if (c.logfile == "-") {
        // use standard out
        openLogFile(stdout);
    } else {
        // use log file
        if (!openLogFile(c.logfile.c_str())) {
            printError("cannot open log file '%s'.", c.logfile.c_str());
            return 1;
        }
    }
    
    setLogLevel(c.verbose);
    
    // print command line options
    if (isLogLevel(LOG_LOW)) {
        printLog(LOG_LOW, "SPIDIR executed with the following arguments:\n");
        for (int i=0; i<argc; i++) {
            printLog(LOG_LOW, "%s ", argv[i]);
        }
        printLog(LOG_LOW, "\n\n");
    }
 
    
    c.printarguments();

    // seed random number generator
    if (c.seed == 0)
        c.seed = time(NULL);
    srand(c.seed);
    printLog(LOG_LOW, "random seed: %d\n", c.seed);
    


    //============================================================
    // read species tree
    SpeciesTree stree_noWGD;

    if (!readNewickTree(c.streefile.c_str(), &stree_noWGD)) {
        printError("error reading species tree '%s'", c.streefile.c_str());
	return 1;
    }

   
    printf("\nSpecies tree from input file, including WGD nodes:\n");
    displayTree(&stree_noWGD, stdout, 0.2);
    writeNewickTree(stdout, &stree_noWGD,true);


    SpeciesTree *WGDstree = removeWGDnodes(&stree_noWGD);


    // write WGD parameters
    printf("\nWGD parameters, if any:\n");
    for (int i=0; i<WGDstree->nWGD; i++) {
       WGDparam * wgd = WGDstree->theWGD[i];
       printf("WGD event number%d:\n", i);
       printf("\tloss probability:  %f\n",wgd->lossProb);
       printf("\tpercent distance:  %f\n",wgd->percentDist);
       printf("\ttotal distance:    %f\n",wgd->totalDist);
       printf("\tnode name at the end of 'before' period: %d\n",
	      wgd->WGD_before->name);
       printf("\tdistance 'before': %f\n", wgd->WGD_before->dist);
       printf("\tnode name at the end of 'at' period:     %d\n",
	      wgd->WGD_at->name);
       printf("\tnode name at the end of 'after' period:  %d\n",
	      wgd->WGD_after->name);
       printf("\tdistance 'after':  %f\n", wgd->WGD_after->dist);
    }

    printf("\nthe WGD species tree\n");
    for (int i=0; i<WGDstree->nnodes; i++) {
      Node *node = WGDstree->nodes[i];
      printf("%d\t%s\n", node->name,(node->longname).c_str());
      for (int i=0; i<node->nchildren; i++){
	printf("\t%d\t%s\n", node->children[i]->name, (node->children[i]->longname).c_str());
      }
    }
    
 

    printf("\nthe WGD species tree graphically:\n");
    writeNewickTree(stdout, WGDstree,true);

    displayTree(WGDstree, stdout, 0.2);  


    printf("\nReduced species tree:\n");
    writeNewickTree(stdout, &stree_noWGD,true);

    displayTree(&stree_noWGD, stdout,0.2);


////////////////////////////////////////////


    stree_noWGD.setDepths();
    
    // read sequences 
    Sequences *aln = readAlignFasta(c.alignfile.c_str());
    auto_ptr<Sequences> aln_ptr(aln);
    if (aln == NULL || !checkSequences(aln->nseqs, aln->seqlen, aln->seqs)) {
        printError("bad alignment file");
        return 1;
    }
    
    // check alignment
    if (aln->nseqs == 0) {
        printError("no sequences");
        return 1;
    }    

    // read SPIDIR parameters
    SpidirParams *params = NULL;
    if (c.paramsfile != "") {
      params = readSpidirParams(c.paramsfile.c_str());
    } else {
        // use default params
        printLog(LOG_LOW, 
         "Note: Rate parameters (-p) were not specified. Using flat prior.\n");
        params = new NullSpidirParams();
    }

    auto_ptr<SpidirParams> params_ptr(params);
    if (params == NULL) {
        printError("error reading parameters file '%s'", c.paramsfile.c_str());
        return 1;
    }
    
    // check params
    if (!params->order(&stree_noWGD)) {
      // Warning: branch 'names' in param file corresponds to branch 'names' in stree_noWGD
      // only if the species tree files with/without WGD use the same left/right.
      // fixit: Place a warning in manual later on.
        printError("parameters do not correspond to the given species tree");
        return 1;
    }        
   
    // determine background base frequency
    float bgfreq[4];

    if (c.bgfreqstr == "") {
        // compute frequency from alignment
        computeBgfreq(aln->nseqs, aln->seqs, bgfreq);
    } else {
        // use supplied frequency
        vector<string> tokens = split(c.bgfreqstr.c_str(), ",");
        if (tokens.size() != 4) {
            printError("bgfreq requires four base frequencies e.g .25,.25,.25,.25");
            return 1;
        }
        for (unsigned int i=0; i<tokens.size(); i++) {
            if (sscanf(tokens[i].c_str(), "%f", &bgfreq[i]) != 1) {
                printError("bgfreq must be floats");
                return 1;
            }
        }
    }
    

    // read gene2species map
    Gene2species mapping;
    if (!mapping.read(c.smapfile.c_str())) {
        printError("error reading gene2species mapping '%s'", 
                   c.smapfile.c_str());
        return 1;
    }

    // get gene names
    ExtendArray<string> genes(0, aln->nseqs);
    genes.extend(aln->names, aln->nseqs);    
    
    // get species names
    ExtendArray<string> species(stree_noWGD.nnodes);
    stree_noWGD.getNames(species);
    
    // make gene to species mapping
    int nnodes = aln->nseqs * 2 - 1; // number of nodes in the plain gene tree

    ExtendArray<int> gene2species(nnodes); 
    mapping.getMap(genes, aln->nseqs, species, stree_noWGD.nnodes, gene2species);
    
    // get initial gene tree by neighbor joining
    Tree *tree = getInitialTree(genes, aln->nseqs, aln->seqlen, aln->seqs,
                                &stree_noWGD, gene2species);
   
    auto_ptr<Tree> tree_ptr(tree);

    printf("\nInitial gene tree from Neighbor-joining:\n");
    fflush(stdout);  
    displayTree(tree);
    // writeNewickTree(stdout, tree,true);     

    //========================================================
    // determine kappa

    if (c.kappa < 0) {
        const float minkappa = .4;
        const float maxkappa = 5.0;
        const float stepkappa = .1;

        printLog(LOG_LOW, "finding optimum kappa...\n");
        // get initial branch lengths
        parsimony(tree, aln->nseqs, aln->seqs); 
        c.kappa =  findMLKappaHky(tree, aln->nseqs, aln->seqs, 
                                  bgfreq, 
                                  minkappa, maxkappa, stepkappa);
        printLog(LOG_LOW, "optimum kappa = %f\n", c.kappa);
    }
    

    //=====================================================
    // init model
    //  Model *model;    


    //me
    SpimapModel *model;

    
    if (WGDstree->nWGD > 0){
      //we will have  some rates for WGD 
      extendRateParamToWGDnodes(params,WGDstree);
    } 

    fflush(stdout);

    SpimapModel *m = new SpimapModel(nnodes, WGDstree, &stree_noWGD, params,
                                     gene2species,
                                     c.pretime, 
                                     c.duprate, 
                                     c.lossrate,
                                     c.priorSamples,
                                     !c.priorExact,
                                     true,c.q);
    /*

     printf("now printing species tree\n");
     fflush(stdout);

     for (int i=0; i<WGDstree->nnodes; i++) {
       Node *node = WGDstree->nodes[i];
       printf("%d\t%s\n", node->name,(node->longname).c_str());
       for (int i=0; i<node->nchildren; i++){
	 printf("\t%d\t%s\n", node->children[i]->name, (node->children[i]->longname).c_str());
       }
     }

     fflush(stdout);

    */

     m->setLikelihoodFunc(new HkySeqLikelihood(
        aln->nseqs, aln->seqlen, aln->seqs, 
        bgfreq, c.kappa, c.lkiter, 
        c.minlen, c.maxlen));

    
    model = m;
    auto_ptr<Model> model_ptr(model);


    //========================================================
    // initialize search
    
    // init topology proposer
    
    float sprrate = .5;


    DefaultSearch prop(c.niter, c.quickiter,
                       &stree_noWGD, gene2species,
                       c.duprate, c.lossrate,
                       sprrate,c.propid);


    MixProposer *proposer = &prop.mix;
    MixProposer *proposer2 = &prop.mix2;

     // init search
    TreeSearchClimb *search = new TreeSearchClimb(model, proposer, proposer2);
 

    auto_ptr<TreeSearch> search_ptr(search);

    // load correct tree
    Tree correctTree;    
    if (c.correctFile != "") {
        if (!readNewickTree(c.correctFile.c_str(), &correctTree)) {
            printError("cannot read correct tree '%s'", c.correctFile.c_str());
            return 1;
        }
        // TODO: aborts if leaves mismatch, should catch error
        correctTree.reorderLeaves(genes);
        proposer->setCorrect(&correctTree);
    }
    
    printf("\nCorrect gene tree from -c file:\n");
    //displayTree(&correctTree);
    //=======================================================
    // search
    time_t startTime = time(NULL);
    // here is when the first reconciliation happens

    Tree *toptree = search->search(tree, genes, 
                                   aln->nseqs, aln->seqlen, aln->seqs, c.outprefix, c.method,c.keepTreeSampled,c.keepDupLoss,c.observingsomething);

    // return 1;

    auto_ptr<Tree> toptree_ptr(toptree);
    if (c.bootiter > 1) {
	if (!bootstrap(aln, genes, search, c.bootiter, c.outprefix))
	    return 1;
    }
    
   
    fflush(stdout);

    //========================================================
    // output final tree
    
    if (isLogLevel(LOG_LOW))
        displayTree(toptree, getLogFile());


    // output recon
    if (c.outputRecon) {
        setInternalNames(toptree);
        setInternalNames(WGDstree);
	string outreconFilename = c.outprefix  + ".recon";	
	  writeRecon(outreconFilename.c_str(), toptree, WGDstree, search->getmodel()->recon, search->getmodel()->events);
    }



    // output gene tree
    string outtreeFilename = c.outprefix  + ".tree";
    writeNewickTree(outtreeFilename.c_str(), toptree);

    

    /////////////////////////////
    printf("now printing species tree\n");
    displayTree(WGDstree, stdout, 0.2);
    fflush(stdout);

    
     for (int i=0; i<WGDstree->nnodes; i++) {
       Node *node = WGDstree->nodes[i];
       printf("%d\t%s\n", node->name,(node->longname).c_str());
       for (int i=0; i<node->nchildren; i++){
	 printf("\t%d\t%s\n", node->children[i]->name, (node->children[i]->longname).c_str());
       }
     }
    
     /////////////////////////////

    

     printf("now printing final gene tree\n");
     displayTree(toptree, stdout,60,2);
  
     
     for (int i=0; i<toptree->nnodes; i++) {
       Node *node = toptree->nodes[i];
       printf("%d\t%s\n", node->name,(node->longname).c_str());
       for (int i=0; i<node->nchildren; i++){
	 printf("\t%d\t%s\n", node->children[i]->name, (node->children[i]->longname).c_str());
       }
     }
     
     
     ///////////////////////////////


     // log tree correctness
     if (c.correctFile != "") {      
       if (proposer->seenCorrect()) {
	 printLog(LOG_LOW, "SEARCH: correct visited\n");
       } else {
	 printLog(LOG_LOW, "SEARCH: correct NEVER visited\n");
       }
        
       if (toptree->sameTopology(&correctTree)) {
	 printLog(LOG_LOW, "RESULT: correct\n");
       } else {
	 printLog(LOG_LOW, "RESULT: wrong\n");
       }
     }
  
    
    fflush(stdout);
    // log runtime
    time_t runtime = time(NULL) - startTime;
    printLog(LOG_LOW, "seq runtime:\t%f\n", model->seq_runtime);
    printLog(LOG_LOW, "branch runtime:\t%f\n", model->branch_runtime);
    printLog(LOG_LOW, "topology runtime:\t%f\n", model->top_runtime);
    printLog(LOG_LOW, "proposal runtime:\t%f\n", search->proposal_runtime);
    printLog(LOG_LOW, "runtime seconds:\t%d\n", runtime);
    printLog(LOG_LOW, "runtime minutes:\t%.1f\n", float(runtime / 60.0));
    printLog(LOG_LOW, "runtime hours:\t%.1f\n", float(runtime / 3600.0));


    
    closeLogFile();
    
}

