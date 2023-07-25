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
  DATA_VECTOR( Kmax );  // max number of occasions
  DATA_MATRIX( log_effort ); //number of trials/surveys per trap. will be repeated by the number of grid cells
  DATA_ARRAY( dist );     //distance from each grid cell center to each trap location
  DATA_ARRAY( y );        //integer matrix of encounter frequencies - rows are individuals, columns are grid cells, and elements are time steps. includes row for all 0 histories
  DATA_INTEGER( surveyarea ); // area of the state space (in km squared)
  DATA_VECTOR( garea_prop ); // proportion of state space covered by each grid cell
  DATA_INTEGER( nstate ); //number of states individuals can be in - needed for 'open' part
  DATA_IVECTOR( year ); // proportion of state space covered by each grid cell
  DATA_IVECTOR( month ); // proportion of state space covered by each grid cell
  DATA_MATRIX( obs_state ); //matrix of known states for tagged individuals - includes 0s for individuals not tagged for indexing purposes
  DATA_IMATRIX( taghist ); //matrix of known states for tagged individuals - includes 0s for individuals not tagged for indexing purposes


  // Parameters
  PARAMETER_VECTOR( alpha0 );      //intercept
  PARAMETER_VECTOR( log_alpha1 );  //constrained to be positive
  PARAMETER( alpha2 );      //coefficient on search effort
  //PARAMETER( alpha3 );      //effect of previous tagging on encounter rates
  PARAMETER_VECTOR( log_lambda );  //intercept term for poisson-distributed N 
  PARAMETER_VECTOR( logit_gamma );  // recruitment - logit scale to constrain from 0-1
  PARAMETER_VECTOR( logit_phi );           // survival
  PARAMETER_VECTOR( logit_imm );           // immigration back into surveyed pop

  //Derived parameters
  vector<Type> alpha1 = exp( log_alpha1 );                //Constrained to be positive
  vector<Type> sigma = pow( Type(1)/(Type(2)*alpha1), Type(0.5) ); //derived sigma 
  vector<Type> lambda = exp( log_lambda ); // number of individuals 
  vector<Type> gamma = exp( logit_gamma )/( 1 + exp( logit_gamma ) );
  vector<Type> phi = exp( logit_phi )/( 1 + exp( logit_phi ) );
  vector<Type> imm = exp( logit_imm )/( 1 + exp( logit_imm ) );

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
  vector<Type> ones_for_T( ntrap ); // for converting logpois to probcap after incorporating previous tagging status
  ones_for_T.fill(1);

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
  array<Type> loginv_probcap(ngrid, ntrap, ntime); // set up array to store 1 minus probability of encounter given trap and individual activity center locations for calculating marginal encounter prob below

  vector<Type> probg( ntrap );
  vector<Type> obs( ngrid );        // vector to store observed encounter histories in likelihood loop
  vector<Type> nocc( ntrap );
  array<Type> posterior( ngrid, nenc, ntime );    // set up matrix to store grib probs surface for each individual 
  matrix<Type> zeroposterior( ngrid, ntime );    // for those with all zero histories

  // Prob of capture in each trap if latent COA in each grid cell in each primary period
  for( int t=0; t<ntime; t++ ){
   for( int g=0; g<ngrid; g++){
     for( int j=0; j<ntrap; j++ ){ 
    //   probcap(g,j,t) = p0 * exp( -alpha1*dist(g,j)*dist(g,j) ); // keep these separate to make clear
      //  probcap(g,j,t) = Type(1) - exp( -exp( alpha0(0) + alpha2*log_effort(j,t) ) * exp( -alpha1(0)*dist(g,j)*dist(g,j) ) ); 
        logpois(g,j,t) = alpha0( 0 ) + alpha2*log_effort(j,t) - alpha1( 0 )*dist(g,j)*dist(g,j);
        probcap(g,j,t) = Type(1) - exp( -exp( logpois(g,j,t) ) );

     // log of 1 minus prob cap for computing marg prob of encounter below  
        loginv_probcap(g,j,t) = log( Type(1) - probcap(g,j,t) ); 
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
      zeroposterior.col(t) = aweight_lcl_zero/marg_lik_zero(t);
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
           posterior.col(t).col(i) = zeroposterior.col(t);
        // if seen
        } else {
         for( int g=0; g<ngrid; g++ ){
           probg = probcap.col(t).matrix().row(g);
           // if including behavioral tagging effect - hash above probg row and unhash next two
           //probg = logpois.col(t).matrix().row(g); 
           //probg = ones_for_T - exp( -exp( probg + alpha3*taghist(i,t)) );               
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// Open (HMM) part of the model
// Set up TPM
// storage
int Tplus1 = ntime + 1;
int Tminus1 = ntime - 1;

// initial distribution of Markov chain
vector<Type> inidist(3);
inidist(0) = 1 - gamma(0);
inidist(1) = gamma(0);
inidist(2) = 0;

//Probability 'born' derived from gamma
// Derive vector of complementary probability sequence
vector<Type> log_cprob( gamma.size() );
Type log_qgamma = Type(0);
// version from BPA book
for (int n=0; n<gamma.size(); n++) {
      log_cprob(n) = log( gamma(n) ) + log_qgamma;
      log_qgamma = log_qgamma + log( Type(1) - gamma(n) );
    }

vector<Type> cprob = exp(log_cprob);
Type psi = cprob.sum();
vector<Type> b = cprob/psi;

// set up matrices for transitions between states
array<Type> tpm( nstate, nstate, Tminus1); // transition probability matrix - needs to have a row for each state in each time step. cols specify probs between states
// second time step will have no immigrants because no one will have been in state 3
tpm(0,0,0) = 1 - gamma( 1 ); //2 ); // (1);
tpm(0,1,0) = gamma( 1 ); //2 ); // (1);
tpm(0,2,0) = 0;
tpm(1,0,0) = 0;
tpm(1,1,0) = phi(0); //if by month (1); //(if time specific or constant (0);
tpm(1,2,0) = 1 - phi(0); //1); //(0);
tpm(2,0,0) = 0;
tpm(2,1,0) = 0; //0; // if third state is absorbing = death/permanent emigration
tpm(2,2,0) = 1; // 1; if third state is absorbing

for (int t=1; t<Tminus1; t++) {
  tpm(0,0,t) = 1 - gamma(t+1); //( 1 ); //month(t) ); 
  tpm(0,1,t) = gamma(t+1); //( 1 ); //month(t) ); 
  tpm(0,2,t) = 0;
  tpm(1,0,t) = 0;
  tpm(1,1,t) = phi( t );  //if constant (0); if time-specific (t); if by month ( month(t)-1 )  
  tpm(1,2,t) = 1 - phi( t ); //(t); //month(t)-1 ); //; 
  tpm(2,0,t) = 0;
  tpm(2,1,t) = imm( t-1 ); //same as for phi but if time-specific (t-1)
  tpm(2,2,t) = 1 - imm( t-1 ); //( 0 );
 }

/// set up storage/temporary objects for HMM state loop
// matrix for initial states - everyone not recruited at t=0
  matrix<Type> gam(ntime, nstate); // 'belief state' for each latent state in each time period
  vector<Type> acc(nstate); // 'accumulator' of evidence for each state based on obs
  matrix<Type> tpmsub(nstate,nstate);
  Type sum_sprob; // storage for sum of likelihoods for forward algorithm
 // matrix<Type> enc_lik(ntrap, nstate);
  vector<Type> enc_lik(nstate); // prob of encounter given state
  vector<Type> cumenc_lik(nstate); // cumulative prob of encounter given state
  vector<Type> na_enc_lik(nstate); //vector for when survey skipped
  na_enc_lik.fill(1); // fill that vector with 1s

  matrix<Type> probf( ngrid, ntrap );
  //matrix<Type> ones_for_M( ngrid, ntrap ); // for converting logpois to probcap after incorporating previous tagging status
  //ones_for_M.fill(1);

// for each individual encountered
for( int i=0; i<nenc; i++) {
  // reset everything to 0
  gam.setZero();

// calculate maginal encounter probability if present given SCR parameter estimates for each time step
for (int t=0; t<ntime; t++) {

  //for( int g=0; g<ngrid; g++){
  //   for( int j=0; j<ntrap; j++ ){ 
           // if including behavioral tagging effect 
  //         probf(g,j) = logpois(g,j,t) + alpha3*taghist(i,t); 
  //         probf(g,j) = Type(1) - exp( -exp( probf(g,j) ) );

     // log of 1 minus prob cap for computing marg prob of encounter below  
   //       loginv_probcap(g,j,t) = log( Type(1) - probf(g,j) );   
   //     }
   //    }

    vector<Type> loginv_t = loginv_probcap.col(t).matrix().rowwise().sum();
    // calc rowsums
    vector<Type> ps_marg = ( 1 - exp( loginv_t ) ) * garea_prop;
    // calc marginal prob of not encountering
    Type margzenc = 1-ps_marg.sum();
    
    //calc transition probs into each state based on state at time t-1
    if( t==0 ){ // would have had to be seen to be tagged in first time step
    // for first time step, tpm has already been applied via use of initial distribution
      // just have to multiply that by prob of encounter informed by obs
       // if not seen
       if( y.col(t).matrix().row(i).sum() == 0 ){
         enc_lik(0) = 1; // not encountered with 100% prob if not recruited
         enc_lik(1) = pow( margzenc, Kmax(t) ); //pow( marg_lik_zero(t), Kmax(t) ); // prob present but not seen
         enc_lik(2) = 1; // not encountered with 100% prob if dead
       // if seen
       } else {
        enc_lik(0) = 0; // can't be seen if not recruited 
        enc_lik(1) = 1 - pow( margzenc, Kmax(t) );// pow( marg_lik_zero(t), Kmax(t)); // present and seen at all with prob specified above
        enc_lik(2) = 0; // can't be seen if dead
        } 
      // multiply prob of state times prob of detection to give prob of encounter
      // using initial distributio and encounter likelihood   
      gam.row(t) = inidist * enc_lik; 

     // for subsequent time steps
     } else { 
      // subset out tpm
     tpmsub = tpm.col(t-1).matrix();
      // if state unknown in previous time step use tpm
     for (int j=0; j<nstate; j++) { // j is current state
        for( int k=0; k<nstate; k++) { // k is previous state
             acc(k) = gam(t-1, k) * tpmsub(k,j); //belief state times transition prob

        if( y.col(t).matrix().row(i).sum() == 0 ){
         enc_lik(0) = 1; // not encountered with 100% prob if not recruited
         enc_lik(1) = pow( margzenc, Kmax(t) ); //pow( marg_lik_zero(t), Kmax(t) ); // prob present but not seen
         enc_lik(2) = 1; // not encountered with 100% prob if dead
       // if seen
       } else {
         enc_lik(0) = 0; // can't be seen if not recruited 
         enc_lik(1) = 1 - pow( margzenc, Kmax(t) );// pow( marg_lik_zero(t), Kmax(t)); // present and seen at all with prob specified above
         enc_lik(2) = 0; // can't be seen if dead
        } 
      // multiply prob of state times prob of detection to give prob of encounter
      acc(k) = acc(k) * enc_lik(j); 
     } 
    gam(t,j) = acc.sum();
    }
   } }
     // contribution to the likelihood
     // note you only need to use last time step here because forward probability implies last iteration retains information of all intermediate state probs
    jnll_comp(0) -= log( gam.row( ntime-1 ).sum() ); 
  } 

  //////////////////////////
  // Derive N based on estimates
  Type Nold = lambda(0)/b(0);
  // To estimate N
  vector<Type> B( ntime ); // recruits from M
  vector<Type> Ntot( ntime ); // total that have traveled through in each time step
  // Pop in first time step is equal to B0
  B(0) = lambda(0);
  Ntot(0) = B(0);
 // Use parameter estimates to update in each step after that
 // Use parameter estimates to update in each step after that
// for (int t=1; t<ntime; t++) { 
 //  if( t==1 ){
 //    B(t) = lambda(t) - lambda(t-1)*phi( t-1 ); // can't immigrate back in because no one is in state 3 yet
 //   } else {
  //   B(t) = lambda(t) - lambda(t-1)*phi( t-1 ) - ( Ntot(t-1) - lambda(t-1) )*imm( t-2 ); 
 //   }
 //   Ntot(t) = Ntot(t-1) + B(t);
 // }
  Type Nsup = Ntot(ntime-1);
  Type Dsup = Nsup/surveyarea;

  /////// Reporting
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

  REPORT( posterior );
  REPORT( zeroposterior );

  // Derived parameters
  ADREPORT( alpha0 );
  ADREPORT( sigma );
  ADREPORT( alpha2 );
  //ADREPORT( alpha3 );
  ADREPORT( b );
  ADREPORT( phi );
  ADREPORT( imm );

  ADREPORT( lambda );
  //ADREPORT( D );
  ADREPORT( Nsup );
  ADREPORT( Nold );
  // ADREPORT( B );
  //ADREPORT( Dsup );

  return jnll;
}

