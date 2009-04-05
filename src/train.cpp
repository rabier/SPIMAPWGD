/*=============================================================================

    SPIDIR    
    Parameter Estimation (Training)
    
    mldist.cpp
    started: Fri Mar 27 14:22:58 EDT 2009

=============================================================================*/

// c++ headers
#include <math.h>
#include <time.h>

// 3rd party
#include <gsl/gsl_multimin.h>

// spidir headers
#include "common.h"
#include "Matrix.h"
#include "seq_likelihood.h"
#include "parsimony.h"
#include "spidir.h"
#include "Tree.h"


namespace spidir {



class RatesEM
{
public:

    RatesEM(int ntrees, int nspecies, int nrates, 
            float **lengths, float *times,
            float *sp_alpha, float *sp_beta, 
            float gene_alpha, float gene_beta) : 
        ntrees(ntrees),
        nspecies(nspecies),          
        lengths(lengths),
        times(times),
        sp_alpha(sp_alpha),
        sp_beta(sp_beta),
        gene_alpha(gene_alpha),
        gene_beta(gene_beta),
        nrates(nrates),
        gtab(ntrees, nrates),
        pgtab(ntrees, nrates)
    {
        // allocate optimizer
        const int ndim = 2;
        opt = gsl_multimin_fdfminimizer_alloc(
            gsl_multimin_fdfminimizer_vector_bfgs2, ndim);

        // setup optimizer for gene rates
        opt_gene_rate.f = &gene_rate_f;
        opt_gene_rate.df = &gene_rate_df;
        opt_gene_rate.fdf = &gene_rate_fdf;
        opt_gene_rate.n = ndim;
        opt_gene_rate.params = this; 

        // setup optimizer for species rates
        opt_sp_rate.f = &sp_rate_f;
        opt_sp_rate.df = &sp_rate_df;
        opt_sp_rate.fdf = &sp_rate_fdf;
        opt_sp_rate.n = ndim;
        opt_sp_rate.params = this;
    }


    ~RatesEM()
    {
        gsl_multimin_fdfminimizer_free(opt);
    }

    // gene rates function and derivative

    static double gene_rate_f(const gsl_vector *x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_alpha = gsl_vector_get(x, 0);
        double gene_beta = gsl_vector_get(x, 1);

        // clamp gamma params
        if (gene_alpha < .001 ||
            gene_beta < .001)
            return -INFINITY;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                float lngamma = gammalog(em->gtab[j][k], gene_alpha, gene_beta);
                if (!isnan(lngamma))
                    sum += em->pgtab[j][k] * lngamma;
            }
        }

        return -sum;
    }
    

    static void gene_rate_df(const gsl_vector *x, void *params, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_alpha = gsl_vector_get(x, 0);
        double gene_beta = gsl_vector_get(x, 1);

        // clamp gamma params
        if (gene_alpha < .001)
            gene_alpha = .001;
        if (gene_beta < .001)
            gene_beta = .001;

        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                float g = em->gtab[j][k];
                float lngamma = gammalog(g, gene_alpha, gene_beta);
                float gamma = exp(lngamma);

                if (!isnan(gamma) && gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] * 
                        gammaDerivA(g, gene_alpha, gene_beta) / gamma;
                    beta_sum += em->pgtab[j][k] * 
                        gammaDerivB(g, gene_alpha, gene_beta) / gamma;
                }
            }
        }


        //printf(". g = (%f, %f)\n", 
        //       gsl_vector_get(x, 0),
        //       gsl_vector_get(x, 1));
        //printf(". dg = (%f, %f)\n", alpha_sum, beta_sum);

        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);
    }

    static void gene_rate_fdf(const gsl_vector *x, void *params, 
                              double *f, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double gene_alpha = gsl_vector_get(x, 0);
        double gene_beta = gsl_vector_get(x, 1);     
   
        
        // clamp gamma params
        if (gene_alpha < .001)
            gene_alpha = .001;
        if (gene_beta < .001)
            gene_beta = .001;        

        //printf("try %f %f\n", gene_alpha, gene_beta);

        double sum = 0.0;
        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                float g = em->gtab[j][k];
                float lngamma = gammalog(g, gene_alpha, gene_beta);
                float gamma = exp(lngamma);
                
                sum += em->pgtab[j][k] * lngamma;
                
                // TODO: substitute a better test
                if (gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] *
                        gammaDerivA(g, gene_alpha, gene_beta) / gamma;
                    beta_sum += em->pgtab[j][k] * 
                        gammaDerivB(g, gene_alpha, gene_beta) / gamma;
                }
            }
        }

        // set return
        *f = -sum;

        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);

        //printf("gene f = %f; fdf = (%f, %f); a=%f, b=%f\n", 
        //       sum,
        //       alpha_sum / em->nrates, 
        //       beta_sum / em->nrates,
        //       gene_alpha, gene_beta);
    }


    //======================================================
    // gene rates function and derivative

    static double sp_rate_f(const gsl_vector *x, void *params)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);
        int i = em->cur_species;

        double sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                float lngamma = gammalog(em->lengths[j][i], sp_alpha_i, 
                                sp_beta_i / (em->gtab[j][k] *
                                             em->times[i]));
                if (!isnan(lngamma))
                    sum += em->pgtab[j][k] * lngamma;
            }
        }

        return -sum;
    }
    

    static void sp_rate_df(const gsl_vector *x, void *params, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);
        int i = em->cur_species;

        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                float g = em->gtab[j][k];
                float bgt = sp_beta_i / (g * em->times[i]);
                float l = em->lengths[j][i];
                float lngamma = gammalog(l, sp_alpha_i, bgt);
                float gamma = exp(lngamma);

                if (gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] * 
                        gammaDerivA(l, sp_alpha_i, bgt) / gamma;
                    beta_sum += em->pgtab[j][k] *
                        gammaDerivB(l, sp_alpha_i, bgt) /
                        (gamma * g * em->times[i]);
                }
            }
        }

        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);
    }

    static void sp_rate_fdf(const gsl_vector *x, void *params, 
                              double *f, gsl_vector *df)
    {
        RatesEM *em = (RatesEM*) params;
        double sp_alpha_i = gsl_vector_get(x, 0);
        double sp_beta_i = gsl_vector_get(x, 1);  
        int i = em->cur_species;
   
        double sum = 0.0;
        double alpha_sum = 0.0;
        double beta_sum = 0.0;
        for (int j=0; j<em->ntrees; j++) {
            for (int k=0; k<em->nrates; k++) {
                float g = em->gtab[j][k];
                float bgt = sp_beta_i / (g * em->times[i]);
                float l = em->lengths[j][i];
                float lngamma = gammalog(l, sp_alpha_i, bgt);
                float gamma = exp(lngamma);

                sum += em->pgtab[j][k] * lngamma;
                
                if (gamma != 0.0) {
                    alpha_sum += em->pgtab[j][k] *
                        gammaDerivA(l, sp_alpha_i, bgt) / gamma;
                    beta_sum += em->pgtab[j][k] * 
                        gammaDerivB(l, sp_alpha_i, bgt) /
                        (gamma * g * em->times[i]);
                }
            }
        }

        // set return
        *f = -sum;

        //printf("sp %d, sum=%f, %f, %f\n", i, sum, alpha_sum, beta_sum);

        gsl_vector_set(df, 0, -alpha_sum);
        gsl_vector_set(df, 1, -beta_sum);
    }


    float likelihood()
    {
        float logl = 0.0;

        EStep();

        for (int j=0; j<ntrees; j++) {
            float sum = 0.0;
            for (int k=0; k<nrates; k++) {
                float prod = 0.0;
                for (int i=0; i<nspecies; i++)
                    prod += gammalog(lengths[j][i], sp_alpha[i], 
                                     sp_beta[i] / (times[i] * gtab[j][k]));
                sum += pgtab[j][k] * exp(prod);
            }
            logl += log(sum);
        }

        return logl;
    }

    //====================

    void init_params()
    {
        for (int i=0; i<nspecies; i++) {
            sp_alpha[i] = 1.0;
            sp_beta[i] = 1.0;
        }

        gene_alpha = 1.0;
        gene_beta = 1.0;
    }

    inline double gene_post(double g, double A, double B, double C)
    {
        return pow(g, B) * exp(A*g - C / g);
    }
    

    float find_upper_g(float m, float A, float B, float C,
                       float tol1=.05, float tol2=.01)
    {        
        double fm = gene_post(m, A, B, C);
        float top = 2*m;
        float bot = m;

        // extent top
        while (true) {
            double ftop = gene_post(top, A, B, C);
            if (ftop/fm <= tol2)
                break;
            top *= 2;
        }

        // binary search
        while (true) {
            float u = (top + bot) / 2.0;
            double fu = gene_post(u, A, B, C);

            if (fu / fm > tol1)
                bot = u;
            else if (fu / fm < tol2)
                top = u;
            else
                return u;
        }
    }


    float find_lower_g(float m, float A, float B, float C,
                       float tol1=.05, float tol2=.01)
    {
        double fm = gene_post(m, A, B, C);
        float top = m;
        float bot = 0;
            
        // binary search
        while (true) {
            float u = (top + bot) / 2.0;
            double fu = gene_post(u, A, B, C);

            if (fu / fm > tol1)
                top = u;
            else if (fu / fm < tol2)
                bot = u;
            else
                return u;
        }
    }


    // populate gene rate posteriors
    void EStep()
    {
        // temp variables for PDF of posterior gene rate
        float x[nrates+1];
        double y[nrates+1];

        // determine commonly used coefficients
        float A = - gene_beta;
        float B = gene_alpha - 1.0;
        for (int i=0; i<nspecies; i++)
            B -= sp_alpha[i];

        for (int j=0; j<ntrees; j++) {

            // determine commonly used coefficients
            float C = 0.0;
            for (int i=0; i<nspecies; i++)
                C += sp_beta[i] * lengths[j][i] / times[i];

            // find main range of gene rates
            float mid = (-B - sqrt(B*B - 4*A*C)) / (2*A);
            float top = find_upper_g(mid, A, B, C);
            float bot = find_lower_g(mid, A, B, C);

            int half_nrates = (nrates + 1) / 2;
            float step1 = (mid - bot) / half_nrates;
            float step2 = (top - mid) / (nrates + 1 - half_nrates);

            // compute x, y for posterior gene rate PDF
            for (int k=0; k<half_nrates; k++) {
                x[k] = bot + step1 * k;
                y[k] = gene_post(x[k], A, B, C);
            }
            for (int k=half_nrates; k<nrates+1; k++) {
                x[k] = mid + step2 * (k-half_nrates);
                y[k] = gene_post(x[k], A, B, C);
            }

            // compute gtab and pgtab
            double total = 0.0;
            for (int k=0; k<nrates; k++) {
                gtab[j][k] = (x[k] + x[k+1]) / 2.0;
                pgtab[j][k] = (y[k] + y[k+1]) * (x[k+1] - x[k]) / 2.0;
                total += pgtab[j][k];
            }

            // normalize pgtab[j]
            for (int k=0; k<nrates; k++) {
                pgtab[j][k] /= total;
            }
        }
    }
    

    // maximize each model parameter given the hidden data estimated from
    // last iteration
    void MStep()
    {
        // optimization config
        double step_size = .1;
        double tol = .1;
        const double epsabs = .001;
        gsl_vector *init_x = gsl_vector_alloc(2);
        int status;
        
        
        // optimize gene rate parameters
        gsl_vector_set(init_x, 0, gene_alpha);
        gsl_vector_set(init_x, 1, gene_beta);
        gsl_multimin_fdfminimizer_set(opt, &opt_gene_rate, init_x, 
                                      step_size, tol);      
        do {
            // do one iteration
            status = gsl_multimin_fdfminimizer_iterate(opt);
            if (status)
                break;        
            // get gradient
            status = gsl_multimin_test_gradient(opt->gradient, epsabs);
            
        } while (status == GSL_CONTINUE);

        gene_alpha = gsl_vector_get(opt->x, 0);
        gene_beta = gsl_vector_get(opt->x, 1);

        //printf("g = (%f, %f)\n", gene_alpha, gene_beta); 
        
        // optimize each species rate parmater set
        for (int i=0; i<nspecies; i++) {
            cur_species = i;
            gsl_vector_set(init_x, 0, sp_alpha[i]);
            gsl_vector_set(init_x, 1, sp_beta[i]);
            gsl_multimin_fdfminimizer_set(opt, &opt_sp_rate, init_x, 
                                          step_size, tol);            
            do {
                // do one iteration
                status = gsl_multimin_fdfminimizer_iterate(opt);
                if (status)
                    break;        
                // get gradient
                status = gsl_multimin_test_gradient(opt->gradient, epsabs);
            } while (status == GSL_CONTINUE);

            sp_alpha[i] = gsl_vector_get(opt->x, 0);
            sp_beta[i] = gsl_vector_get(opt->x, 1);

            //printf("sp[%d] = (%f, %f)\n", i, sp_alpha[i], sp_beta[i]);
        }

        gsl_vector_free(init_x);
    }

    // data
    int ntrees;
    int nspecies;
    float **lengths;

    // given fixed parameters
    float *times;

    // model parameters
    float *sp_alpha;
    float *sp_beta;
    float gene_alpha;
    float gene_beta;
    
    //protected:

    // hidden data
    int nrates;
    Matrix<float> gtab;
    Matrix<double> pgtab;

    // optimizer
    gsl_multimin_fdfminimizer *opt;
    gsl_multimin_function_fdf opt_gene_rate;
    gsl_multimin_function_fdf opt_sp_rate;
    int cur_species;
};


extern "C" {

void train(int ntrees, int nspecies, float **lengths, float *times,
           float *sp_alpha, float *sp_beta, float *gene_alpha, float *gene_beta,
           int nrates, int max_iter)
{

    /*
    for (int i=0; i<ntrees; i++) {
        printFloatArray(lengths[i], nspecies);
        printf("\n");
    }
    
    printf("times ");
    printFloatArray(times, nspecies);
    printf("\n");
    */

    RatesEM em(ntrees, nspecies, nrates, lengths, times,
               sp_alpha, sp_beta, *gene_alpha, *gene_beta);
    
    
    // make initial guess for model parameters
    em.init_params();
    
    // iterate until convergence
    for (int iter=0; iter<max_iter; iter++) {
        em.EStep();
        em.MStep();

        //printf("logl: %f\n", em.likelihood());
    }
      
}


RatesEM *allocRatesEM(int ntrees, int nspecies, int nrates,
                      float **lengths, float *times,
                      float *sp_alpha, float *sp_beta, 
                      float gene_alpha, float gene_beta)
{
    // copy lengths
    float **lengths2 = new float* [ntrees];
    for (int j=0; j<ntrees; j++) {
        lengths2[j] = new float [nspecies];
        for (int i=0; i<nspecies; i++)
            lengths2[j][i] = lengths[j][i];
    }

    float *times2 = new float [nspecies];
    float *sp_alpha2 = new float [nspecies];
    float *sp_beta2 = new float [nspecies];

    for (int i=0; i<nspecies; i++) {
        times2[i] = times[i];
        sp_alpha2[i] = sp_alpha[i];
        sp_beta2[i] = sp_beta[i];
    }

    return new RatesEM(ntrees, nspecies, nrates, lengths2, times2,
                              sp_alpha2, sp_beta2, gene_alpha, gene_beta);
}


void freeRatesEM(RatesEM *em)
{
    for (int j=0; j<em->ntrees; j++)
        delete [] em->lengths[j];
    delete [] em->lengths;

    delete [] em->times;
    delete [] em->sp_alpha;
    delete [] em->sp_beta;
    
    delete em;
}


void RatesEM_Init(RatesEM *em)
{
    em->init_params();
}


void RatesEM_EStep(RatesEM* em)
{
    em->EStep();
}

void RatesEM_MStep(RatesEM* em)
{
    em->MStep();
}

float RatesEM_likelihood(RatesEM *em)
{
    /*
    for (int j=0; j<em->ntrees; j++) {
        for (int k=0; k<em->nrates; k++) {
            printf("%f ", em->gtab[j][k]);
        }
        printf("\n");
        }*/

    return em->likelihood();
}


void RatesEM_getParams(RatesEM *em, float *params)
{
    params[0] = em->gene_alpha;
    params[1] = em->gene_beta;

    for (int i=0; i<em->nspecies; i++) {
        params[2+2*i] = em->sp_alpha[i];
        params[2+2*i+1] = em->sp_beta[i];
    }
}



} // extern "C"



} // spidir
