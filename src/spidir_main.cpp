/*=============================================================================

    Test SPIDIR functions

=============================================================================*/

// c++ headers
#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <vector>

// third party
#include <gsl/gsl_errno.h>

// spidir headers
#include "common.h"
#include "phylogeny.h"
#include "parsimony.h"
#include "search.h"
#include "Matrix.h"
#include "ConfigParam.h"
#include "Sequences.h"
#include "spidir.h"


#define VERSION_INFO  "\
   ___    SPIDIR v2.0 2009 \n\
  /0 0\\   SPecies Informed DIstanced-base Reconstruction \n\
  \\___/   Matt Rasmussen \n\
 /// \\\\\\  CSAIL, MIT \n\
"


using namespace std;
using namespace spidir;

class SpidirConfig
{
public:

    SpidirConfig() 
    {}

    int parseArgs(int argc, char **argv)
    {
	// parse arguments
	ConfigParser config;
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
		   ("-p", "--param", "<spidir params file>", &paramsfile, 
		    "SPIDIR branch length parameters file"));
	config.add(new ConfigParam<string>
		   ("-o", "--output", "<output filename prefix>", 
		    &outprefix, "spidir",
		    "prefix for all output filenames"));
    
    
	config.add(new ConfigParamComment("Sequence evolution model"));
	config.add(new ConfigParam<float>
		   ("-k", "--kappa", "<transition/transversion ratio>", 
		    &kappa, 1.0,
		    "used for HKY model (default=1.0)"));
	config.add(new ConfigParam<string>
		   ("-f", "--bgfreq", "<A freq>,<C ferq>,<G freq>,<T freq>", 
		    &bgfreqstr, "",
		    "background frequencies (default=0.25,0.25,0.25,0.25)"));

	config.add(new ConfigParamComment("Dup/loss evolution model"));
	config.add(new ConfigParam<string>
		   ("", "--prior", "hky|spidir|duploss|birthdeath|none", 
		    &prioropt, "spidir",
		    "function for prior (default=spidir)"));
	config.add(new ConfigParam<float>
		   ("-D", "--duprate", "<duplication rate>", 
		    &duprate, 0.1,
		    "rate of a gene duplication (default=0.1)"));
	config.add(new ConfigParam<float>
		   ("-L", "--lossrate", "<loss rate>", 
		    &lossrate, 0.1,
		    "probability of loss (default=0.1)"));
	config.add(new ConfigParam<float>
		   ("-P", "--pretime", "<pre-speciation time parameter>", 
		    &pretime, 1.00,
		    "lambda param of pre-speciation distribution (default=1.0)"));

    
	config.add(new ConfigParamComment("Miscellaneous"));
	config.add(new ConfigParam<string>
		   ("", "--search", "climb|mcmc", 
		    &search, "climb", 
		    "search algorithm (default=climb)"));
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
	config.add(new ConfigParam<string>
		   ("-c", "--correct", "<correct tree file>", &correctFile, ""
		    "check if correct tree is visited in search"));
	config.add(new ConfigParam<int>
		   ("", "--prior_samples", "<number of samples>",
		    &priorSamples, 100,
		    "number of samples to use in branch prior integration (default: 100)"));
	config.add(new ConfigSwitch
		   ("", "--prior_exact", 
		    &priorExact,
		    "Use an exact calculation of branch prior"));
	config.add(new ConfigSwitch
		   ("-g", "--gene_rate", &estGenerate, "estimate generate"));
        config.add(new ConfigParam<int>
                   ("", "--quicksamples", "<number of quick samples>",
                    &quickSamples, 10,
                    "number of samples in quick search"));
	config.add(new ConfigSwitch
		   ("", "--no_spr_nbr", &noSprNbr, "do not use SPR NBR"));


	config.add(new ConfigParamComment("Information"));
	config.add(new ConfigParam<int>
		   ("-V", "--verbose", "<verbosity level>", 
		    &verbose, LOG_LOW, 
		    "verbosity level 0=quiet, 1=low, 2=medium, 3=high"));
	config.add(new ConfigParam<string>
		   ("", "--log", "<log filename>", &logfile, "", 
		    "log filename.  Use '-' to display on stdout."));
	config.add(new ConfigSwitch
		   ("-v", "--version", &version, "display version information"));
	config.add(new ConfigSwitch
		   ("-h", "--help", &help, 
		    "display help information"));

    
    
	if (!config.parse(argc, (const char**) argv)) {
	    if (argc < 2)
		config.printHelp();
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

    string alignfile;    
    string smapfile;
    string streefile;
    string paramsfile;
    string outprefix;
    string search;
    string correctFile;
    string prioropt;
    int niter;
    float kappa;
    string bgfreqstr;
    float pretime;
    float duprate;
    float lossrate;
    int priorSamples;
    bool priorExact;
    string logfile;
    int verbose;
    bool help;
    bool version;
    bool estGenerate;
    int bootiter;
    int quickiter;
    int quickSamples;
    bool noSprNbr;
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

	for (int i=1; i<bootiter; i++) {
	    printLog(LOG_LOW, "bootstrap %d of %d\n", i, bootiter);
	    resampleAlign(aln, &aln2);
            
	    boottree = search->search(NULL, genes, 
				      aln2.nseqs, aln2.seqlen, aln2.seqs);

	    boottree->setLeafNames(genes);
	    boottree->writeNewick(bootfile, NULL, 0, true);
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
    // seed random number generator
    srand(time(NULL));

    SpidirConfig c;
    int ret = c.parseArgs(argc, argv);
    if (ret)
	return ret;
    
    //============================================================
    // output filenames
    string outtreeFilename = c.outprefix  + ".tree";
    
    // use default log filename
    if (c.logfile == "")
        c.logfile = c.outprefix + ".log";
    
    if (c.logfile == "-") {
        // use standard out
        openLogFile(stdout);
    } else {
        if (!openLogFile(c.logfile.c_str())) {
            printError("cannot open log file '%s'.", c.logfile.c_str());
            return 1;
        }
    }
    
    setLogLevel(c.verbose);
    
    if (isLogLevel(LOG_LOW)) {
        printLog(LOG_LOW, "SPIDIR executed with the following arguments:\n");
        for (int i=0; i<argc; i++) {
            printLog(LOG_LOW, "%s ", argv[i]);
        }
        printLog(LOG_LOW, "\n\n");
    }
    
    
    //============================================================
    // read species tree
    SpeciesTree stree;
    if (!stree.readNewick(c.streefile.c_str())) {
        printError("error reading species tree '%s'", c.streefile.c_str());
        return 1;
    }
    stree.setDepths();
    
    
    // read sequences
    Sequences *aln;
    
    if ((aln = readAlignFasta(c.alignfile.c_str())) == NULL ||
        !checkSequences(aln->nseqs, aln->seqlen, aln->seqs)) {
        printError("bad alignment file");
        return 1;
    }
    
    if (aln->nseqs < 3) {
        printError("too few sequences");
        return 1;
    }
    

    // read SPIDIR parameters
    SpidirParams *params;
    if ((params = readSpidirParams(c.paramsfile.c_str())) == NULL)
    {
        printError("error reading parameters file '%s'", c.paramsfile.c_str());
        return 1;
    }
    
    if (!params->order(&stree)) {
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
    
    
    int nnodes = aln->nseqs * 2 - 1;

    // read gene2species map
    Gene2species mapping;
    if (!mapping.read(c.smapfile.c_str())) {
        printError("error reading gene2species mapping '%s'", c.smapfile.c_str());
        return 1;
    }

    // produce mapping array
    ExtendArray<string> genes(0, aln->nseqs);
    genes.extend(aln->names, aln->nseqs);    
    
    ExtendArray<string> species(stree.nnodes);
    stree.getLeafNames(species);
    
    ExtendArray<int> gene2species(nnodes);
    mapping.getMap(genes, aln->nseqs, species, stree.nnodes, gene2species);
    
    
    //=====================================================
    // init prior function
    Prior *prior;    

    if (c.prioropt == "none")
        prior = new Prior();
    
    else if (c.prioropt == "spidir")
        prior = new SpidirPrior(nnodes, &stree, params, 
				gene2species,
				c.pretime, 
				c.duprate, 
				c.lossrate,
				c.priorSamples,
				!c.priorExact,
				c.estGenerate);
    /*
    else if (c.lkfuncopt == "duploss")
        lkfunc = new SpidirBranchLikelihoodFunc(nnodes, &stree, params, 
                                                gene2species,
                                                c.pretime, 
						c.duprate,
						c.lossrate, 
                                                c.estGenerate,
                                                true);
    else if (c.lkfuncopt == "hky") 
        lkfunc = new HkyBranchLikelihoodFunc(aln->nseqs, aln->seqlen, aln->seqs, 
                                             bgfreq, c.tsvratio);
    else if (c.lkfuncopt == "birthdeath") 
        lkfunc = new BranchLikelihoodFunc();
    */
    else {
        printError("unknown prior '%s'", c.prioropt.c_str());
        return 1;
    }
    

    //========================================================
    // branch lengths
    
    // determine branch length algorithm
    BranchLengthFitter *fitter = NULL;
    const int maxiter = 10;
    fitter = new HkyFitter(aln->nseqs, aln->seqlen, aln->seqs, 
			   bgfreq, c.kappa, maxiter);
    
    
    //========================================================
    // initialize search
    
    // init topology proposer
    const int radius = 3;
    TopologyProposer *proposer2;
    //if (c.noSprNbr)
        proposer2 = new SprNniProposer(&stree, gene2species, c.niter);
    //else
    //    proposer2 = new SprNbrProposer(&stree, gene2species, c.niter, radius);

    DupLossProposer proposer(proposer2, &stree, gene2species, 
			     c.duprate, c.lossrate,
                             c.quickiter, c.niter);//, 
    //                             c.quickSamples);
    

    TreeSearch *search = NULL;
    if (c.search == "climb") {
	search = new TreeSearchClimb(prior, &proposer, fitter);
    } else {
        printError("unknown search '%s'", c.search.c_str());
        return 1;
    }

    // load correct tree
    Tree correctTree;    
    if (c.correctFile != "") {
        if (!correctTree.readNewick(c.correctFile.c_str())) {
            printError("cannot read correct tree '%s'", c.correctFile.c_str());
            return 1;
        }
        // TODO: aborts if leaves mismatch, should catch error
        correctTree.reorderLeaves(genes);
        proposer.setCorrect(&correctTree);
    }
       
    
    time_t startTime = time(NULL);
    
    //=======================================================
    // setup gsl
    gsl_set_error_handler_off();


    //=======================================================
    // search
    Tree *toptree;
    
    Tree *tree = getInitialTree(genes, aln->nseqs, aln->seqlen, aln->seqs,
                                &stree, gene2species);

    if (c.bootiter <= 1) {
	toptree = search->search(tree, genes, 
				 aln->nseqs, aln->seqlen, aln->seqs);
    } else {        
	toptree = search->search(tree, genes, 
				 aln->nseqs, aln->seqlen, aln->seqs);
	if (!bootstrap(aln, genes, search, c.bootiter, c.outprefix)) {
	    return 1;
	}
    }
    
    //========================================================
    // output final tree
    
    displayTree(toptree);

    toptree->setLeafNames(genes);
    toptree->writeNewick(outtreeFilename.c_str());
    
    
    // log tree correctness
    if (c.correctFile != "") {
        if (proposer.seenCorrect()) {
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

        
    // log runtime
    time_t runtime = time(NULL) - startTime;
    printLog(LOG_LOW, "seq runtime:\t%f\n", fitter->runtime);
    printLog(LOG_LOW, "branch runtime:\t%f\n", prior->branch_runtime);
    printLog(LOG_LOW, "topology runtime:\t%f\n", prior->top_runtime);
    printLog(LOG_LOW, "runtime seconds:\t%d\n", runtime);
    printLog(LOG_LOW, "runtime minutes:\t%.1f\n", float(runtime / 60.0));
    printLog(LOG_LOW, "runtime hours:\t%.1f\n", float(runtime / 3600.0));
    closeLogFile();
    
    // clean up
    delete toptree;
    delete params;
    delete fitter;
    delete prior;
    delete search;
    delete proposer2;
}


//=============================================================================
// OLD CODE

    /*
    if (c.search == "mcmc") {
        string mcmcfilename = c.outprefix + ".mcmc";
        FILE *mcmcfile = NULL;
        
        if (!(mcmcfile = fopen(mcmcfilename.c_str(), "w"))) {
            printError("cannot open mcmc file '%s'", mcmcfilename.c_str());
            return 1;
        }
        SampleFunc samples(mcmcfile);
        
        if (c.bootiter != 1) {
            printError("Cannot use bootstrap with MCMC");
            return 1;
        }
    
        toptree = searchMCMC(tree, 
                             genes, aln->nseqs, aln->seqlen, aln->seqs,
                             &samples,
                             lkfunc,
                             &proposer,
                             fitter);

    } else 
    */
