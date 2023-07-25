// Fit SCR model to spatial encounter histories
#include <TMB.hpp> // Include TMB-specific macros and functions, including dependencies such as CppAD and Eigen

// Function for detecting NAs
template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

// Start definition of objective function
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace R_inla; // Namespaces are the equivalent of packages in R - includes a bunch of pre-specified functions
  using namespace Eigen;
  using namespace density;

 // Data
  DATA_INTEGER( nenc );   //number of individuals ever encountered
  DATA_IVECTOR( nsurv );  // number of individuals encountered in each survey period 
  DATA_INTEGER( ntrap );  //number of traps (for us, number of grid cells visited during the course of each survey)
  DATA_INTEGER( ngrid );  //number of grid cells
  DATA_INTEGER( ntime );  //number of time steps/trials
  DATA_MATRIX( K );       //number of trials per trap (cols) and time periods (rows)
  DATA_MATRIX( log_effort ); //number of trials/surveys per trap. will be repeated by the number of grid cells
  DATA_ARRAY( dist );     //distance from each grid cell center to each trap location
  DATA_ARRAY( y );        //integer matrix of encounter frequencies - rows are individuals, columns are grid cells, and elements are time steps. includes row for all 0 histories
  DATA_INTEGER( surveyarea ); // area of the state space (in km squared)
  DATA_VECTOR( garea_prop ); // proportion of state space covered by each grid cell
  DATA_IVECTOR( year ); // proportion of state space covered by each grid cell
  DATA_IVECTOR( month ); // proportion of state space covered by each grid cell

  // Parameters
  PARAMETER_VECTOR( alpha0 );      //intercept
  PARAMETER_VECTOR( log_alpha1 );  //constrained to be positive
  PARAMETER( alpha2 );      //coefficient on search effort
  PARAMETER_VECTOR( log_lambda );  //intercept term for poisson-distributed N 

  //Derived parameters
  //Type p0 = exp( alpha0 )/( 1 + exp(alpha0) );    //logit transform for alpha0 to constrain prob between 0 and 1
  //Type p0 = exp( alpha0 );                          //if modeling as number of encounters instead
  vector<Type> alpha1 = exp( log_alpha1 );                //Constrained to be positive
  vector<Type> sigma = pow( Type(1)/(Type(2)*alpha1), Type(0.5) ); //derived sigma 
  vector<Type> lambda = exp( log_lambda ); // number of individuals 

  //Create storage for different components of the likelihood and other bits needed for model fitting
  vector<Type> jnll_comp( 3 ); //storage for negloglik summed over all individual encounter histories (including 0s)
  // includes two components - the combinatorial term and the sum of the marginal likelihoods over all encounter histories
  jnll_comp.setZero();
  // for individual contributions to the nLL
  matrix<Type> loglik_obs( nenc, ntime ); // for encountered aka observed individuals
  loglik_obs.setZero();
  vector<Type> loglik_zero( ntime ); // for zeros
  loglik_zero.setZero();
  
  // objects needed for likelihood loops
  matrix<Type> y_zero( ntime, ntrap ); // for zero encounter histories
  y_zero.setZero(); //set to zero
  vector<Type> y_zero_T( ntrap ); // for subsetting in loop

  // conditional and marginal likelihood intermediates (need these because TMB loses dim of an array when taking the log or exp)
  matrix<Type> log_cond_lik( nenc, ngrid );
  vector<Type> log_cond_lik_zero( ngrid );
  matrix<Type> cond_lik( nenc, ngrid );
  vector<Type> tmp( ntrap ); // storage for loglik in each grid cell for lik loop
  vector<Type> lcl( ngrid );
  vector<Type> lcl_zero( ngrid );
  vector<Type> aweight_lcl( ngrid );
  vector<Type> aweight_lcl_zero( ngrid );
  matrix<Type> marg_lik( nenc, ntime );
  vector<Type> marg_lik_zero( ntime );
  
  // storage of various compononets for model fitting
  array<Type> probcap( ngrid, ntrap, ntime );     // set up array to store probability of encounter given trap and individual activity center locations
  array<Type> logpois( ngrid, ntrap, ntime );
  vector<Type> probg( ntrap );
  vector<Type> obs( ngrid );        // vector to store observed encounter histories in likelihood loop
  vector<Type> nocc( ntrap );
  array<Type> posterior( ngrid, nenc, ntime );    // set up matrix to store grib probs surface for each individual   

  // Prob of capture in each trap if latent COA in each grid cell in each primary period
  for( int t=0; t<ntime; t++ ){
   for( int g=0; g<ngrid; g++){
     for( int j=0; j<ntrap; j++ ){ 
    //   probcap(g,j,t) = p0 * exp( -alpha1*dist(g,j)*dist(g,j) ); // keep these separate to make clear
      //  probcap(g,j,t) = Type(1) - exp( -exp( alpha0(0) + alpha2*log_effort(j,t) ) * exp( -alpha1(0)*dist(g,j)*dist(g,j) ) ); 
        logpois(g,j,t) = alpha0( 0 ) + alpha2*log_effort(j,t) - alpha1(0)*dist(g,j)*dist(g,j);
        probcap(g,j,t) = Type(1) - exp( -exp( logpois(g,j,t) ) );
        }
      }
    }  

  // Encounter history probability for observed individuals if present
   // All zero encounter history probability
   for( int t=0; t<ntime; t++ ){
        nocc = K.row(t);
        y_zero_T = y_zero.row(t);
   for( int g=0; g<ngrid; g++ ){
        probg = probcap.col(t).matrix().row(g);
        tmp = dbinom( y_zero_T, nocc, probg, TRUE );
        log_cond_lik_zero(g) = tmp.sum();
    }
      lcl_zero = log_cond_lik_zero; 
      aweight_lcl_zero = exp(lcl_zero) * garea_prop; //* ( Type(1)/ngrid );
      marg_lik_zero(t) = aweight_lcl_zero.sum(); 
      loglik_zero(t) = log( marg_lik_zero(t) );
  }

  // Non-zero encounter history
  for( int t=0; t<ntime; t++ ){
        nocc = K.row(t);
    for( int i=0; i<nenc; i++ ){
       // subset out obs for individual in that time step
        obs = y.col(t).matrix().row(i);
        // if not seen, don't contribute to the likelihood here
        if( obs.sum() == 0 ) {
           //loglik_obs(i,t) = 0;
           // but can still calc prob of COA for mapping if in area
           posterior.col(t).col(i) = aweight_lcl_zero/marg_lik_zero(t);
        // if seen
        } else {
         for( int g=0; g<ngrid; g++ ){
           probg = probcap.col(t).matrix().row(g);         
           tmp = dbinom( obs, nocc, probg, TRUE );
           log_cond_lik(i,g) = tmp.sum();
          }
        lcl = log_cond_lik.row(i); 
        aweight_lcl = exp(lcl) * garea_prop; //* ( Type(1)/ngrid );
        marg_lik(i,t) = aweight_lcl.sum(); 
        loglik_obs(i,t) = log( marg_lik(i,t) );
        posterior.col(t).col(i) = aweight_lcl/marg_lik(i,t);
        }
       }  
      }
 
 // Poisson integrated likelihood (form used by SECR package- see SCR book Chapter 6)
 for( int t=0; t<ntime; t++ ){
  jnll_comp(0) -= nsurv(t) * log( lambda(t) );
  jnll_comp(1) -= loglik_obs.col(t).sum();
  jnll_comp(2) -= -lambda(t) * ( Type(1) - exp( loglik_zero(t) ) );
}

  // Density in each time step
  vector<Type> D = lambda/surveyarea;

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );

  // checking to make sure things work as intended when developing
  REPORT( obs );
  //REPORT( Pm );
  //REPORT( gridProb );
  REPORT( probcap );
//  REPORT( probg );
  //REPORT( probcap_vec );
  REPORT( aweight_lcl );
  REPORT( loglik_obs );
  REPORT( loglik_zero );
  REPORT( tmp );
  //REPORT( tmp.sum() );
  REPORT( lcl );
  REPORT( marg_lik );
  REPORT( loglik_obs );
  REPORT( marg_lik_zero );

  // Derived parameters
  ADREPORT( alpha0 );
  ADREPORT( sigma );
  ADREPORT( alpha2 );
  ADREPORT( lambda );
  ADREPORT( D );

  return jnll;
}

