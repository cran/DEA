// 
// esta version usa las matrices con la notacion
// habitual en estadistica para las matrices de datos.
//
//   28/01/2008  versión inicial
//   08/02/2008  corregido bug ("no" en lugar de "nD") en modelos dea.ccr.io.mul,
//               dea.ccr.oo.mul, dea.bcc.io.mul, dea.bcc.oo.mul, dea.bcc.oo.env,
//               dea.add.mul
//   10/02/2008  añadida la opción del presolver

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "glpk.h"

/********************************************************/
/*                                                      */
/*  comienza el modelo ccr input orientado en la forma  */
/*  multiplicador                                       */
/*                                                      */
/********************************************************/

void ccr_io_mul (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {

  // comienza la funcion ccr

  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

  // numero de restricciones; debera ser nD+1
  int nR;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s, t;


  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.
  double obje, restri;

  int ne_max =( (ni + no) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no filas y nD+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;
  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar

  for ( i = 1 ; i <= nD ; i++ ) {

    for ( j = 1 ; j <= ni ; j++ ) {

      restri = *( xx + ((j-1) * nD ) + (i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( j = ni+1 ; j <= ni+no ; j++){

      restri = *( yy + ((j-ni-1) * nD ) + (i-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }
  }
  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD filas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.

    /*********     BUCLE GORDO      *********/


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

  // creamos el problema, le asignamos como nombre "dea_ccr_io_mul"
  // e indicamos que es de maximizacion

  lp = lpx_create_prob();
  lpx_set_prob_name(lp, "dea_ccr_io_mul");
  lpx_set_obj_dir(lp,LPX_MAX);

  // añadimos filas, es decir, restricciones, al problema.
  // de momento solo añado las que no cambian

  lpx_add_rows(lp,nD);

  // establecemos el nombre y los limites de cada fila,
  // es decir, de cada variable auxiliar.

  for (t = 1 ; t <= nD ; t++ ) {

    lpx_set_row_name (lp, t, *(nombre_DMUs + (t-1)));
    lpx_set_row_bnds (lp, t, LPX_LO, 0.0, 1000.0);
  }
  
  // esto corresponde a la restriccion de igualdad
  // que la pongo la ultima


  // añadimos columnas, es decir, variables estructurales

  lpx_add_cols ( lp, ni + no);

  // establecemos el nombre y los limites de cada variable estructural

  for ( t = 1 ; t <= ni ; t++ ) {

    lpx_set_col_name (lp, t, *(nombre_inputs + (t-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }

  for (t = ni+1 ; t <= ni+no ; t++ ){

    lpx_set_col_name (lp, t, *(nombre_outputs + (t-ni-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }


    // añado la restriccion de igualdad; nR indica la posicion en la que
    // la añado, que debera ser nD+1

    nR = lpx_add_rows (lp, 1);
    lpx_set_row_name (lp, nR, "normalizacion input virtual");
    lpx_set_row_bnds (lp, nR, LPX_FX, 1.0, 1.0);


    for ( j = 1 ; j <= ni ; j++) {

      obje = 0.0;
      lpx_set_obj_coef (lp, j, obje);

      restri = *( xx + ((j-1) * nD ) + (i-1));

      // rellenamos la ultima restriccion, la nR-esima

      if (restri != 0.0) {
	nf = nf+1; 
	ia[nf] = nR;
	ja[nf] = j;
	ar[nf] = restri;
      }
    }

    for ( j = ni+1 ; j <= ni+no ; j++ ) {

      obje = *(yy + ((j-ni-1)*nD) + (i-1));
      lpx_set_obj_coef (lp, j, obje);
    }


    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex (lp);

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion

    for( j = 1 ; j <= ni+no ; j++ ) {
      *(zz +(i-1) + (j-1)*nD ) = lpx_get_col_prim (lp , j );
    }


    *(zz+ (i-1) + (ni+no)*nD ) = lpx_get_obj_val (lp);

    lpx_delete_prob (lp);

  }


} 

/********************************************************/
/*************** fin de la funcion ccr_io_mul ***********/
/********************************************************/

/********************************************************/
/*                                                      */
/*  comienza el modelo ccr input orientado en la forma  */
/*  envolvente                                          */
/*                                                      */
/********************************************************/


void ccr_io_env (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double teta;
  double limi;

  int ne_max =( (ni + no) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no filas y nD+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }
  }

  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD columnas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

    // creamos el problema, le asignamos como nombre "dea_ccr_io_env"
    // e indicamos que es de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name(lp, "dea_ccr_io_env");
    lpx_set_obj_dir(lp,LPX_MIN);

    // añadimos filas, es decir, restricciones, al problema.

    lpx_add_rows(lp, ni + no);
    
    // establecemos el nombre y los limites de cada fila,
    // es decir, de cada variable auxiliar.

    for (t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp, t, *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp, t, LPX_LO, 0.0, 1000.0);
    }

    for (t = ni+1 ; t <= ni + no ; t++ ) {

      lpx_set_row_name (lp, t, *(nombre_outputs + (t-ni-1)));
      lpx_set_row_bnds (lp, t, LPX_LO, *( yy + (i-1) + (t-ni-1)*nD ), 1000.0);
    }


    // añadimos columnas, es decir, variables estructurales
    // van a ser nD+1

    lpx_add_cols ( lp, nD+1);


    // establecemos el nombre y los limites de cada variable estructural

    for ( t = 1 ; t <= nD ; t++ ) {

      lpx_set_col_name (lp, t, strcat( ss , *(nombre_DMUs + (t-1)) ) );
      strcpy (ss , tt );
      lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

    }


    // meto la ultima columna; la llamo "teta"

    lpx_set_col_name(lp, nD+1, "teta");
    lpx_set_col_bnds(lp, nD+1, LPX_FR, 0.0, 0.0);

    // relleno esa columna

    for(j = 1; j <= ni ; j++ ) {

      restri = *(xx + (i-1) + (j-1)*nD);
      if (restri != 0.0 ) {
	nf = nf +1 ;
	ia[nf] = j ;
	ja[nf] = nD + 1;
	ar[nf] = restri ;
      }
    }



    // creo la funcion objetivo

    for (j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef(lp , j , 0.0 );
    }

    lpx_set_obj_coef(lp, nD+1 , 1.0);


    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    // elimino mensajes coñazo
    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex ( lp );

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d; phase I: OK, %c %c %c %c; ", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d; phase I: FAULT; ", i );
    }; 

    // extraigo la solucion

    teta = lpx_get_obj_val (lp);
    *(zz + (i-1) + nD * nD ) = teta ;

    lpx_delete_prob (lp);

    /*********  comienza la segunda etapa ****************/

    nf = ne;

    // creo el problema, le asigno nombre e indico que es
    // de maximixacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_ccr_io_env_2");
    lpx_set_obj_dir (lp , LPX_MAX);

    // añado las filas, que van a ser ni+no, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no );

    for ( t = 1 ; t <= ni ; t++ ) {
      limi = *(xx + (i-1) + (t-1) * nD);
      limi = (-1) * teta * limi;
      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {
      limi = *(yy + (i-1) + (t-ni-1) * nD );
      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    //completo la matriz de restricciones

    nf = ne ;

    for (j = 1 ; j <= ni+no ; j++ ) {
      nf = nf+1;
      ia[nf] = j ;
      ja[nf] = j + nD ;
      ar[nf] = -1.0 ;
    }

    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 1.0 );
    }

    // cargo la matriz de restricciones

    lpx_load_matrix (lp , nf , ia , ja , ar );

    // elimino mensajes coñazo

    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex ( lp );

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("phase II: OK, %c %c %c %c \n", kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("phase II: FAULT \n" );
    }; 

    // extraigo la solucion
    // solo me interesan el valor de la funcion objetivo y
    // las lambdas. El valor de la funcion objetivo se va 
    // a la ultima columna de zz

    *(zz + (i-1) + nD*(nD+1) ) = lpx_get_obj_val(lp);

    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j);
    }

    // borro el problema


    lpx_delete_prob(lp);
    
    /***********  acaba la segunda etapa   ***************/


  }
  // fin del bucle gordo

} 

/*********************************************************/
/**************  fin de la funcion ccr_io_env ************/
/*********************************************************/


/*********************************************************/
/*                                                       */
/*  comienza el modelo ccr output orientado en la forma  */
/*  multiplicador                                        */
/*                                                       */
/*********************************************************/

void ccr_oo_mul (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {

  // comienza la funcion ccr

  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

  // numero de restricciones; debera ser nD+1
  int nR;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s, t;


  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;

  int ne_max =( (ni + no) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no filas y nD+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar

  for ( i = 1 ; i <= nD ; i++ ) {

    for ( j = 1 ; j <= ni ; j++ ) {

      restri = *( xx + ((j-1) * nD ) + (i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( j = ni+1 ; j <= ni+no ; j++){

      restri = *( yy + ((j-ni-1) * nD ) + (i-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }
  }
  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD filas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.

    /*********     BUCLE GORDO      *********/


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

  // creamos el problema, le asignamos como nombre "dea_ccr_oo_mul"
  // e indicamos que es de maximizacion

  lp = lpx_create_prob();
  lpx_set_prob_name(lp, "dea_ccr_oo_mul");
  lpx_set_obj_dir(lp,LPX_MIN);

  // añadimos filas, es decir, restricciones, al problema.
  // (de momento solo añado las que no cambian); establezco
  // tambien el nombre y los limites de cada fila

  lpx_add_rows(lp,nD);

  for (t = 1 ; t <= nD ; t++ ) {

    lpx_set_row_name (lp, t, *(nombre_DMUs + (t-1)));
    lpx_set_row_bnds (lp, t, LPX_LO, 0.0, 1000.0);
  }
  


  // añadimos columnas, es decir, variables estructurales
  // establezco tambien el nombre y los limites de cada una

  lpx_add_cols ( lp, ni + no);

  for ( t = 1 ; t <= ni ; t++ ) {

    lpx_set_col_name (lp, t, *(nombre_inputs + (t-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }

  for (t = ni+1 ; t <= ni+no ; t++ ){

    lpx_set_col_name (lp, t, *(nombre_outputs + (t-ni-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }


    // añado la restriccion de igualdad; nR indica la posicion en la que
    // la añado, que debera ser nD+1

    nR = lpx_add_rows (lp, 1);
    lpx_set_row_name (lp, nR, "normalizacion output virtual");
    lpx_set_row_bnds (lp, nR, LPX_FX, 1.0, 1.0);

    // relleno la restriccion de igualdad y la funcion objetivo

    for ( j = 1 ; j <= ni ; j++) {

      obje = *( xx + (i-1) + (j-1)*nD );
      lpx_set_obj_coef (lp, j, obje);
    }
    

    for ( j = ni+1 ; j <= ni+no ; j++ ) {

      obje = 0.0 ;
      lpx_set_obj_coef (lp, j, obje);

      restri = *( yy + (i-1) + ((j-ni-1) * nD ) );

      if (restri != 0.0) {
	nf = nf+1; 
	ia[nf] = nR;
	ja[nf] = j;
	ar[nf] = restri;
      }
    }


    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    // elimino mensajes coñazo
    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 0 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV );

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // flag del scaling
    //lpx_set_int_parm ( lp , LPX_K_SCALE , 2 );


    // simplex iterations limit
    //lpx_set_int_parm ( lp , LPX_K_ITLIM , 200 );


    // aplico el simplex

    rr = lpx_simplex (lp);

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion

    for( j = 1 ; j <= ni+no ; j++ ) {
      *(zz +(i-1) + (j-1)*nD ) = lpx_get_col_prim (lp , j );
    }


    *(zz+ (i-1) + (ni+no)*nD ) = lpx_get_obj_val (lp);

    lpx_delete_prob (lp);

  }

} 

/*********************************************************/
/*********** fin de la funcion ccr_oo_mul ****************/
/*********************************************************/


/*********************************************************/
/*                                                       */
/*  comienza el modelo ccr output orientado en la forma  */
/*  envolvente                                           */
/*                                                       */
/*********************************************************/


void ccr_oo_env (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double teta;
  double limi;

  int ne_max =( (ni + no) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no filas y nD+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala; 

  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }
  }

  
  // esta parte de la matriz de restricciones (la que corresponde
  // a las primeras nD columnas) ya no va a cambiar. El numero de
  // elementos no nulos que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

    // creamos el problema, le asignamos como nombre "dea_ccr_oo_env"
    // e indicamos que es de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name(lp, "dea_ccr_oo_env");
    lpx_set_obj_dir(lp,LPX_MAX);

    // añadimos filas, es decir, restricciones, al problema.

    lpx_add_rows(lp, ni + no);
    
    // establecemos el nombre y los limites de cada fila,
    // es decir, de cada variable auxiliar.

    for (t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp, t, *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp, t, LPX_LO, -(*(xx + (i-1) + (t-1)*nD)), 1000.0);
    }

    for (t = ni+1 ; t <= ni + no ; t++ ) {

      lpx_set_row_name (lp, t , *(nombre_outputs + (t-ni-1)));
      lpx_set_row_bnds (lp, t , LPX_LO , 0.0 , 1000.0);
    }


    // añadimos columnas, es decir, variables estructurales
    // van a ser nD+1

    lpx_add_cols ( lp, nD+1);


    // establecemos el nombre y los limites de cada variable estructural

    for ( t = 1 ; t <= nD ; t++ ) {

      lpx_set_col_name (lp, t, strcat( ss , *(nombre_DMUs + (t-1)) ) );
      strcpy (ss , tt );
      lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

    }


    // meto la ultima columna; la llamo "teta"
    // y la relleno en la matriz de restricciones

    lpx_set_col_name(lp, nD+1, "teta");
    lpx_set_col_bnds(lp, nD+1, LPX_FR, 0.0, 0.0);

    for(j = ni+1 ; j <= ni+no ; j++ ) {

      restri = *(yy + (i-1) + (j-ni-1)*nD);
      if (restri != 0.0 ) {
	nf = nf +1 ;
	ia[nf] = j ;
	ja[nf] = nD + 1;
	ar[nf] = -restri ;
      }
    }



    // creo la funcion objetivo

    for (j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef(lp , j , 0.0 );
    }

    lpx_set_obj_coef(lp, nD+1 , 1.0);


    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    // elimino mensajes coñazo
    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex ( lp );
    
    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d; phase I: OK, %c %c %c %c; ", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d; phase I: FAULT; ", i );
    }; 

    // extraigo la solucion

    teta = lpx_get_obj_val (lp);
    *(zz + (i-1) + nD * nD ) = teta ;

    lpx_delete_prob (lp);

    /*********  comienza la segunda etapa ****************/

    nf = ne;

    // creo el problema, le asigno nombre e indico que es
    // de maximizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_ccr_oo_env_2");
    lpx_set_obj_dir (lp , LPX_MAX);

    // añado las filas, que van a ser ni+no, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no );

    for ( t = 1 ; t <= ni ; t++ ) {
      limi = *(xx + (i-1) + (t-1) * nD);
      limi = (-1) * limi;
      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {
      limi = *(yy + (i-1) + (t-ni-1) * nD );
      limi = teta*limi;
      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    //completo la matriz de restricciones

    nf = ne ;

    for (j = 1 ; j <= ni+no ; j++ ) {
      nf = nf+1;
      ia[nf] = j ;
      ja[nf] = j + nD ;
      ar[nf] = -1.0 ;
    }

    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 1.0 );
    }

    // cargo la matriz de restricciones

    lpx_load_matrix (lp , nf , ia , ja , ar );

    // elimino mensajes coñazo

    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex ( lp );

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf (" phase II: OK, %c %c %c %c \n" , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf (" phase II: FAULT \n" );
    }; 

    // extraigo la solucion
    // solo me interesan el valor de la funcion objetivo y
    // las lambdas. El valor de la funcion objetivo se va 
    // a la ultima columna de zz

    *(zz + (i-1) + nD*(nD+1) ) = lpx_get_obj_val(lp);

    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j);
    }

    // borro el problema


    lpx_delete_prob(lp);
    
    /***********  acaba la segunda etapa   ***************/


  }
  // fin del bucle gordo

} 

/*********************************************************/
/**************  fin de la funcion ccr_oo_env ************/
/*********************************************************/

/*********************************************************/
/*                                                       */
/*  comienza el modelo bcc input orientado en la forma   */
/*  multiplicador                                        */
/*                                                       */
/*********************************************************/

void bcc_io_mul (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,  
		 int    * information,
		 double * zz) {



  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

  // numero de restricciones; debera ser nD+1
  int nR;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s, t;


  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.
  double obje, restri;

  int ne_max =( (nD+1) * (ni+no+1) );
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene nD+1 filas y ni+no+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala; 

  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar

  for ( i = 1 ; i <= nD ; i++ ) {

    for ( j = 1 ; j <= ni ; j++ ) {

      restri = *( xx + ((j-1) * nD ) + (i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( j = ni+1 ; j <= ni+no ; j++){

      restri = *( yy + ((j-ni-1) * nD ) + (i-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }

    ne=ne+1;
    ia[ne]=i;
    ja[ne]=ni+no+1;
    ar[ne]=1;
  }
  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD filas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.

    /*********     BUCLE GORDO      *********/


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

  // creamos el problema, le asignamos como nombre "dea_bcc_io_mul"
  // e indicamos que es de maximizacion

  lp = lpx_create_prob();
  lpx_set_prob_name(lp, "dea_bcc_io_mul");
  lpx_set_obj_dir(lp,LPX_MAX);

  // añadimos filas, es decir, restricciones, al problema.
  // de momento solo añado las que no cambian

  lpx_add_rows(lp,nD);

  // establecemos el nombre y los limites de cada fila,
  // es decir, de cada variable auxiliar.

  for (t = 1 ; t <= nD ; t++ ) {

    lpx_set_row_name (lp, t, *(nombre_DMUs + (t-1)));
    lpx_set_row_bnds (lp, t, LPX_LO, 0.0, 1000.0);
  }
  
  
  // añadimos columnas, es decir, variables estructurales
  // y establezco los nombres y limites de cada una de ellas

  lpx_add_cols ( lp, ni + no + 1);


  for ( t = 1 ; t <= ni ; t++ ) {

    lpx_set_col_name (lp, t, *(nombre_inputs + (t-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }

  for (t = ni+1 ; t <= ni+no ; t++ ){

    lpx_set_col_name (lp, t, *(nombre_outputs + (t-ni-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }

  lpx_set_col_name (lp , ni+no+1 , "uCero" );
  lpx_set_col_bnds (lp , ni+no+1 , LPX_FR , 0.0 , 0.0 );


    // añado la restriccion de igualdad; nR indica la posicion en la que
    // la añado, que debera ser nD+1

    nR = lpx_add_rows (lp, 1);
    lpx_set_row_name (lp, nR, "normalizacion input virtual");
    lpx_set_row_bnds (lp, nR, LPX_FX, 1.0, 1.0);


    for ( j = 1 ; j <= ni ; j++) {

      obje = 0.0;
      lpx_set_obj_coef (lp, j, obje);

      restri = *( xx + ((j-1) * nD ) + (i-1));

      // rellenamos la ultima restriccion, la nR-esima

      if (restri != 0.0) {
	nf = nf+1; 
	ia[nf] = nR;
	ja[nf] = j;
	ar[nf] = restri;
      }
    }

    for ( j = ni+1 ; j <= ni+no ; j++ ) {

      obje = *(yy + ((j-ni-1)*nD) + (i-1));
      lpx_set_obj_coef (lp, j, obje);
    }

    lpx_set_obj_coef ( lp , ni+no+1 , -1.0 );

    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex (lp);

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion

    for( j = 1 ; j <= ni+no ; j++ ) {
      *(zz +(i-1) + (j-1)*nD ) = lpx_get_col_prim (lp , j );
    }


    *(zz+ (i-1) + (ni+no)*nD ) = lpx_get_obj_val (lp);

    lpx_delete_prob (lp);

  }


} 

/********************************************************/
/*************** fin de la funcion bcc_io_mul ***********/
/********************************************************/


/********************************************************/
/*                                                      */
/*  comienza el modelo bcc input orientado en la forma  */
/*  envolvente                                          */
/*                                                      */
/********************************************************/


void bcc_io_env (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double teta;
  double limi;

  int ne_max =( (nD+1) * (ni+no+1) );
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene nD+1 filas y ni+no+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    ne = ne+1;
    ia[ne] = ni+no+1;
    ja[ne] = j;
    ar[ne] = 1;
  }

  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD columnas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

    // creamos el problema, le asignamos como nombre "dea_bcc_io_env"
    // e indicamos que es de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name(lp, "dea_bcc_io_env");
    lpx_set_obj_dir(lp,LPX_MIN);

    // añadimos filas, es decir, restricciones, al problema.

    lpx_add_rows(lp, ni + no + 1);
    
    // establecemos el nombre y los limites de cada fila,
    // es decir, de cada variable auxiliar.

    for (t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp, t, *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp, t, LPX_LO, 0.0, 1000.0);
    }

    for (t = ni+1 ; t <= ni + no ; t++ ) {

      lpx_set_row_name (lp, t, *(nombre_outputs + (t-ni-1)));
      lpx_set_row_bnds (lp, t, LPX_LO, *( yy + (i-1) + (t-ni-1)*nD ), 1000.0);
    }

    lpx_set_row_name (lp , ni+no+1 , "igual_a_uno" );
    lpx_set_row_bnds (lp , ni+no+1 , LPX_FX , 1.0 , 1.0 );


    // añadimos columnas, es decir, variables estructurales
    // van a ser nD+1

    lpx_add_cols ( lp, nD+1);


    // establecemos el nombre y los limites de cada variable estructural

    for ( t = 1 ; t <= nD ; t++ ) {

      lpx_set_col_name (lp, t, strcat( ss , *(nombre_DMUs + (t-1)) ) );
      strcpy (ss , tt );
      lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

    }


    // meto la ultima columna; la llamo "teta"

    lpx_set_col_name(lp, nD+1, "teta");
    lpx_set_col_bnds(lp, nD+1, LPX_FR, 0.0, 0.0);

    // relleno esa columna

    for(j = 1; j <= ni ; j++ ) {

      restri = *(xx + (i-1) + (j-1)*nD);
      if (restri != 0.0 ) {
	nf = nf +1 ;
	ia[nf] = j ;
	ja[nf] = nD + 1;
	ar[nf] = restri ;
      }
    }



    // creo la funcion objetivo

    for (j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef(lp , j , 0.0 );
    }

    lpx_set_obj_coef(lp, nD+1 , 1.0);


    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    // elimino mensajes coñazo
    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex ( lp );

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d; phase I: OK, %c %c %c %c; ", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d; phase I: FAULT; ", i );
    }; 

    // extraigo la solucion

    teta = lpx_get_obj_val (lp);
    *(zz + (i-1) + nD * nD ) = teta ;

    lpx_delete_prob (lp);

    /*********  comienza la segunda etapa ****************/

    nf = ne;

    // creo el problema, le asigno nombre e indico que es
    // de maximizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_bcc_io_env_2");
    lpx_set_obj_dir (lp , LPX_MAX);

    // añado las filas, que van a ser ni+no+1, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no+1 );

    for ( t = 1 ; t <= ni ; t++ ) {
      limi = *(xx + (i-1) + (t-1) * nD);
      limi = (-1) * teta * limi;
      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {
      limi = *(yy + (i-1) + (t-ni-1) * nD );
      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    lpx_set_row_name (lp , ni+no+1 , "igual_a_uno" );
    lpx_set_row_bnds (lp , ni+no+1 , LPX_FX , 1.0 , 1.0 );

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    //completo la matriz de restricciones

    nf = ne ;

    for (j = 1 ; j <= ni+no ; j++ ) {
      nf = nf+1;
      ia[nf] = j ;
      ja[nf] = j + nD ;
      ar[nf] = -1.0 ;
    }

    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 1.0 );
    }

    // cargo la matriz de restricciones

    lpx_load_matrix (lp , nf , ia , ja , ar );

    // elimino mensajes coñazo

    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex ( lp );

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("phase II: OK, %c %c %c %c \n" , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("phase II: FAULT \n" );
    }; 

    // extraigo la solucion
    // solo me interesan el valor de la funcion objetivo y
    // las lambdas. El valor de la funcion objetivo se va 
    // a la ultima columna de zz

    *(zz + (i-1) + nD*(nD+1) ) = lpx_get_obj_val(lp);

    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j);
    }

    // borro el problema


    lpx_delete_prob(lp);
    
    /***********  acaba la segunda etapa   ***************/


  }
  // fin del bucle gordo

} 


/********************************************************/
/*************** fin de la funcion bcc_io_env ***********/
/********************************************************/


/*********************************************************/
/*                                                       */
/*  comienza el modelo bcc output orientado en la forma  */
/*  multiplicador                                        */
/*                                                       */
/*********************************************************/

void bcc_oo_mul (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;


  // numero de restricciones; debera ser nD+1
  int nR;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s, t;


  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;

  int ne_max =( (nD+1) * (ni + no + 1) );
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene nD+1 filas y ni+no+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar

  for ( i = 1 ; i <= nD ; i++ ) {

    for ( j = 1 ; j <= ni ; j++ ) {

      restri = *( xx + ((j-1) * nD ) + (i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( j = ni+1 ; j <= ni+no ; j++){

      restri = *( yy + ((j-ni-1) * nD ) + (i-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }
    ne = ne+1;
    ia[ne] = i;
    ja[ne] = ni+no+1;
    ar[ne] = -1.0;
  }
  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD filas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.

    /*********     BUCLE GORDO      *********/


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

  // creamos el problema, le asignamos como nombre "dea_bcc_oo_mul"
  // e indicamos que es de maximizacion

  lp = lpx_create_prob();
  lpx_set_prob_name(lp, "dea_bcc_oo_mul");
  lpx_set_obj_dir(lp,LPX_MIN);

  // añadimos filas, es decir, restricciones, al problema.
  // (de momento solo añado las que no cambian); establezco
  // tambien el nombre y los limites de cada fila

  lpx_add_rows(lp,nD);

  for (t = 1 ; t <= nD ; t++ ) {

    lpx_set_row_name (lp, t, *(nombre_DMUs + (t-1)));
    lpx_set_row_bnds (lp, t, LPX_LO, 0.0, 1000.0);
  }
  


  // añadimos columnas, es decir, variables estructurales
  // establezco tambien el nombre y los limites de cada una

  lpx_add_cols ( lp, ni + no + 1);

  for ( t = 1 ; t <= ni ; t++ ) {

    lpx_set_col_name (lp, t, *(nombre_inputs + (t-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }

  for (t = ni+1 ; t <= ni+no ; t++ ){

    lpx_set_col_name (lp, t, *(nombre_outputs + (t-ni-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

  }

  lpx_set_col_name (lp , ni+no+1 , "vCero" );
  lpx_set_col_bnds (lp , ni+no+1 , LPX_FR , 0.0 , 0.0 );


    // añado la restriccion de igualdad; nR indica la posicion en la que
    // la añado, que debera ser nD+1

    nR = lpx_add_rows (lp, 1);
    lpx_set_row_name (lp, nR, "normalizacion output virtual");
    lpx_set_row_bnds (lp, nR, LPX_FX, 1.0, 1.0);

    // relleno la restriccion de igualdad y la funcion objetivo

    for ( j = 1 ; j <= ni ; j++) {

      obje = *( xx + (i-1) + (j-1)*nD );
      lpx_set_obj_coef (lp, j, obje);
    }
    

    for ( j = ni+1 ; j <= ni+no ; j++ ) {

      obje = 0.0 ;
      lpx_set_obj_coef (lp, j, obje);

      restri = *( yy + (i-1) + ((j-ni-1) * nD ) );

      if (restri != 0.0) {
	nf = nf+1; 
	ia[nf] = nR;
	ja[nf] = j;
	ar[nf] = restri;
      }
    }

    lpx_set_obj_coef (lp , ni+no+1 , -1.0 );

    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    // elimino mensajes coñazo
    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);
    
    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // pinto el problema para chequearlo
    // if ( i == 1 )  lpx_print_prob ( lp , "/home/teresita/pru1.txt" ); 

    // aplico el simplex
    rr = lpx_simplex ( lp );
    
    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

     /*
     if ( lpx_get_status (lp)  == LPX_OPT )
       Rprintf (" OK, \n" );
     else
       Rprintf (" FAULT, \n" );
     */

    // extraigo la solucion

    for( j = 1 ; j <= ni+no ; j++ ) {
      *(zz +(i-1) + (j-1)*nD ) = lpx_get_col_prim (lp , j );
    }


    *(zz+ (i-1) + (ni+no)*nD ) = lpx_get_obj_val (lp);

    lpx_delete_prob (lp);

  }

} 

/*********************************************************/
/*********** fin de la funcion bcc_oo_mul ****************/
/*********************************************************/


/*********************************************************/
/*                                                       */
/*  comienza el modelo bcc output orientado en la forma  */
/*  envolvente                                           */
/*                                                       */
/*********************************************************/


void bcc_oo_env (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,	 
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s, t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double teta;
  double limi;

  int ne_max =( (ni + no + 1) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no+1 filas y nD+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 

  //int n_it = 20;

 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    ne = ne + 1;
    ia[ne] = ni+no+1;
    ja[ne] = j;
    ar[ne] = 1;
  }

  
  // esta parte de la matriz de restricciones (la que corresponde
  // a las primeras nD columnas) ya no va a cambiar. El numero de
  // elementos no nulos que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

    // creamos el problema, le asignamos como nombre "dea_bcc_oo_env"
    // e indicamos que es de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name(lp, "dea_bcc_oo_env");
    lpx_set_obj_dir(lp,LPX_MAX);


    // añadimos filas, es decir, restricciones, al problema.

    lpx_add_rows(lp, ni + no + 1);
    
    // establecemos el nombre y los limites de cada fila,
    // es decir, de cada variable auxiliar.

    for (t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp, t, *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp, t, LPX_LO, -(*(xx + (i-1) + (t-1)*nD)), 1000.0);
    }

    for (t = ni+1 ; t <= ni + no ; t++ ) {

      lpx_set_row_name (lp, t , *(nombre_outputs + (t-ni-1)));
      lpx_set_row_bnds (lp, t , LPX_LO , 0.0 , 1000.0);
    }

    lpx_set_row_name (lp , ni+no+1 , "igual_a_uno" );
    lpx_set_row_bnds (lp , ni+no+1 , LPX_FX , 1.0 , 1.0 );
    // añadimos columnas, es decir, variables estructurales
    // van a ser nD+1

    lpx_add_cols ( lp, nD+1);


    // establecemos el nombre y los limites de cada variable estructural

    for ( t = 1 ; t <= nD ; t++ ) {

      lpx_set_col_name (lp, t, strcat( ss , *(nombre_DMUs + (t-1)) ) );
      strcpy (ss , tt );
      lpx_set_col_bnds (lp, t, LPX_LO, 0.0, 1000.0);

    }


    // meto la ultima columna; la llamo "teta"
    // y la relleno en la matriz de restricciones

    lpx_set_col_name(lp, nD+1, "teta");
    lpx_set_col_bnds(lp, nD+1, LPX_FR, 0.0, 0.0);

    for(j = ni+1 ; j <= ni+no ; j++ ) {

      restri = *(yy + (i-1) + (j-ni-1)*nD);
      if (restri != 0.0 ) {
	nf = nf +1 ;
	ia[nf] = j ;
	ja[nf] = nD + 1;
	ar[nf] = -restri ;
      }
    }



    // creo la funcion objetivo

    for (j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef(lp , j , 0.0 );
    }

    lpx_set_obj_coef(lp, nD+1 , 1.0);


    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    lpx_reset_parms (lp );

    // elimino mensajes coñazo
    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );
    
    // pinto el problema para chequearlo
    //if ( i == 1 )  lpx_print_prob ( lp , "/home/teresita/pru1.txt" ); 

    // establezco número máximo de iteraciones
    //lpx_set_int_parm ( lp , LPX_K_ITLIM , n_it );

    // aplico el simplex
    rr = lpx_simplex ( lp );

   if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d; phase I: OK, %c %c %c %c; ", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d; phase II: FAULT; ", i );
    }; 

    // chequeo
    //Rprintf("primera etapa; %d \n",i);
    

    // extraigo la solucion

    teta = lpx_get_obj_val (lp);
    *(zz + (i-1) + nD * nD ) = teta ;

    // chequeo teta
    //Rprintf("teta: %f \n", teta);

    lpx_delete_prob (lp);

    /*********  comienza la segunda etapa ****************/

    nf = ne;

    // creo el problema, le asigno nombre e indico que es
    // de maximizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_bcc_oo_env_2");
    lpx_set_obj_dir (lp , LPX_MAX);



    // añado las filas, que van a ser ni+no+1, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no+1 );

    for ( t = 1 ; t <= ni ; t++ ) {
      limi = *(xx + (i-1) + (t-1) * nD);
      limi = (-1) * limi;
      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {
      limi = *(yy + (i-1) + (t-ni-1) * nD );
      limi = teta*limi;
      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    lpx_set_row_name (lp , ni+no+1 , "igual_a_uno" );
    lpx_set_row_bnds (lp , ni+no+1 , LPX_FX , 1.0 , 1.0 );

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    //completo la matriz de restricciones

    nf = ne ;

    for (j = 1 ; j <= ni+no ; j++ ) {
      nf = nf+1;
      ia[nf] = j ;
      ja[nf] = j + nD ;
      ar[nf] = -1.0 ;
    }

    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 1.0 );
    }

    // cargo la matriz de restricciones

    lpx_load_matrix (lp , nf , ia , ja , ar );

    lpx_reset_parms (lp );

    // elimino mensajes coñazo
    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // chequeo
    //lpx_print_prob ( lp , "/home/teresita/pru*.txt" );


    // otro chequeo
    //Rprintf("segunda etapa; %d \n",i);

    // establezco número máximo de iteraciones
    //lpx_set_int_parm ( lp , LPX_K_ITLIM , n_it );

    // aplico el simplex
    rr = lpx_simplex ( lp );

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("phase II: OK, %c %c %c %c \n", kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("phase II: FAULT \n" );
    }; 

    //Rprintf("ok tras simplex\n");

    // extraigo la solucion
    // solo me interesan el valor de la funcion objetivo y
    // las lambdas. El valor de la funcion objetivo se va 
    // a la ultima columna de zz

    *(zz + (i-1) + nD*(nD+1) ) = lpx_get_obj_val(lp);

    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j);
    }

    // borro el problema


    lpx_delete_prob(lp);
    
    /***********  acaba la segunda etapa   ***************/


  }
  // fin del bucle gordo

} 

/*********************************************************/
/**************  fin de la funcion bcc_oo_env ************/
/*********************************************************/



/*********************************************************/
/*                                                       */
/* comienza el modelo aditivo en la forma multiplicador  */
/*                                                       */
/*                                                       */
/*********************************************************/

void add_mul (char  ** nombre_DMUs,
	      int    * numero_DMUs,
	      char  ** nombre_inputs,
	      int    * numero_inputs,
	      char  ** nombre_outputs,
	      int    * numero_outputs,
	      double * xx,
	      double * yy,
	      int    * presolver,
	      int    * dual_simplex,
	      int    * information,
	      double * zz) {



  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s, t;


  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.
  double obje, restri;

  int ne_max =( (nD) * (ni+no+1) );
  // contiene el numero maximo de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene nD filas y ni+no+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar

  for ( i = 1 ; i <= nD ; i++ ) {

    for ( j = 1 ; j <= ni ; j++ ) {

      restri = *( xx + ((j-1) * nD ) + (i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( j = ni+1 ; j <= ni+no ; j++){

      restri = *( yy + ((j-ni-1) * nD ) + (i-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }

    ne=ne+1;
    ia[ne]=i;
    ja[ne]=ni+no+1;
    ar[ne]=1;
  }
  
  // esta  matriz de restricciones (que corresponde a nD filas
  // y ni+no+1 columnas) ya no va a cambiar. El numero de elementos
  // no nulos que tiene es ne.

    /*********     BUCLE GORDO      *********/


  for ( i = 1 ; i <= nD ; i++ ){

    nf = ne;

  // creamos el problema, le asignamos como nombre "dea_add_mul"
  // e indicamos que es de minimizacion

  lp = lpx_create_prob();
  lpx_set_prob_name(lp, "dea_add_mul");
  lpx_set_obj_dir(lp,LPX_MIN);

  // añadimos las nDfilas, es decir, restricciones, al problema
  // y establecemos el nombre y los limites de cada fila,
  // es decir, de cada variable auxiliar.

  lpx_add_rows(lp,nD);

  for (t = 1 ; t <= nD ; t++ ) {

    lpx_set_row_name (lp, t, *(nombre_DMUs + (t-1)));
    lpx_set_row_bnds (lp, t, LPX_LO, 0.0, 1000.0);
  }
  
  
  // añadimos columnas, es decir, variables estructurales
  // y establezco los nombres y limites de cada una de ellas

  lpx_add_cols ( lp, ni + no + 1);


  for ( t = 1 ; t <= ni ; t++ ) {

    lpx_set_col_name (lp, t, *(nombre_inputs + (t-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 1.0, 1000.0);

  }

  for (t = ni+1 ; t <= ni+no ; t++ ){

    lpx_set_col_name (lp, t, *(nombre_outputs + (t-ni-1)));
    lpx_set_col_bnds (lp, t, LPX_LO, 1.0, 1000.0);

  }

  lpx_set_col_name (lp , ni+no+1 , "uCero" );
  lpx_set_col_bnds (lp , ni+no+1 , LPX_FR , 0.0 , 0.0 );


  //establezco la funcion objetivo

    for ( j = 1 ; j <= ni ; j++) {

      obje = *(xx + (i-1) + (j-1) * nD);
      lpx_set_obj_coef (lp, j, obje);

    }

    for ( j = ni+1 ; j <= ni+no ; j++ ) {

      obje = *(yy + ((j-ni-1)*nD) + (i-1));
      lpx_set_obj_coef (lp, j, -obje);
    }

    lpx_set_obj_coef ( lp , ni+no+1 , 1.0 );

    // cargo ahora la matriz de restricciones
    lpx_load_matrix(lp, nf, ia, ja, ar);

    lpx_set_int_parm ( lp , LPX_K_MSGLEV , 1);

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex (lp);

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion

    for( j = 1 ; j <= ni+no ; j++ ) {
      *(zz +(i-1) + (j-1)*nD ) = lpx_get_col_prim (lp , j );
    }


    *(zz+ (i-1) + (ni+no)*nD ) = lpx_get_obj_val (lp);

    lpx_delete_prob (lp);

  }   // fin del bucle gordo


} 

/********************************************************/
/*************** fin de la funcion add_mul **************/
/********************************************************/


/********************************************************/
/*                                                      */
/*    comienza el modelo add en la forma envolvente     */
/*                                                      */
/*                                                      */
/********************************************************/


void add_env (char  ** nombre_DMUs,
	      int    * numero_DMUs,
	      char  ** nombre_inputs,
	      int    * numero_inputs,
	      char  ** nombre_outputs,
	      int    * numero_outputs,
	      double * xx,
	      double * yy,
	      int    * presolver,
	      int    * dual_simplex,
	      int    * information,
	      double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double limi;

  int ne_max =( (nD+1) * (ni+no+1) );
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene nD+1 filas y ni+no+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = -restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    ne = ne+1;
    ia[ne] = ni+no+1;
    ja[ne] = j;
    ar[ne] = 1;
  }

  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD columnas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  for ( i = 1 ; i <= nD ; i++ ){


    nf = ne;

    // creo el problema, le asigno nombre e indico que es
    // de maximizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "add_env");
    lpx_set_obj_dir (lp , LPX_MAX);

    // añado las filas, que van a ser ni+no+1, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no+1 );

    for ( t = 1 ; t <= ni ; t++ ) {
      limi = *(xx + (i-1) + (t-1) * nD);
      limi = (-1) * limi;
      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {
      limi = *(yy + (i-1) + (t-ni-1) * nD );
      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , limi , 1000.0 );
    }

    lpx_set_row_name (lp , ni+no+1 , "igual_a_uno" );
    lpx_set_row_bnds (lp , ni+no+1 , LPX_FX , 1.0 , 1.0 );

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    //completo la matriz de restricciones

    nf = ne ;

    for (j = 1 ; j <= ni+no ; j++ ) {
      nf = nf+1;
      ia[nf] = j ;
      ja[nf] = j + nD ;
      ar[nf] = -1.0 ;
    }

    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 1.0 );
    }

    // cargo la matriz de restricciones

    lpx_load_matrix (lp , nf , ia , ja , ar );

    // elimino mensajes coñazo

    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex

    rr = lpx_simplex(lp);

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion
    // solo me interesan el valor de la funcion objetivo y
    // las lambdas. El valor de la funcion objetivo se va 
    // a la ultima columna de zz

    *(zz + (i-1) + nD*(nD) ) = lpx_get_obj_val(lp);

    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j);
    }

    // borro el problema


    lpx_delete_prob(lp);


  }
  // fin del bucle gordo

} 


/********************************************************/
/*************** fin de la funcion add_env **************/
/********************************************************/



/********************************************************/
/*                                                      */
/*       comienza el modelo slack based-measure         */
/*                          tipo CCR                    */
/*                                                      */
/********************************************************/


void sbm_ccr (char  ** nombre_DMUs,
	      int    * numero_DMUs,
	      char  ** nombre_inputs,
	      int    * numero_inputs,
	      char  ** nombre_outputs,
	      int    * numero_outputs,
	      double * xx,
	      double * yy,
	      int    * presolver,
	      int    * dual_simplex,
	      int    * information,
	      double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double te;
  double limi;

  int ne_max =( (ni + no) * (nD+3));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no filas y nD+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala; 

  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }
  }


    for (j = 1 ; j <= ni ; j++ ) {
      ne = ne+1;
      ia[ne] = j ;
      ja[ne] = j + nD ;
      ar[ne] = 1.0 ;
    }

    for (j = ni +  1 ; j <= ni+no ; j++ ) {
      ne = ne+1;
      ia[ne] = j ;
      ja[ne] = j + nD ;
      ar[ne] = -1.0 ;
    }


  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD columnas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){


    // creo el problema, le asigno nombre e indico que es
    // de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_sbm");
    lpx_set_obj_dir (lp , LPX_MIN);

    // añado las filas, que van a ser ni+no+1, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no+1 );

    for ( t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , 0.0 , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , 0.0 , 1000.0 );
    }

    lpx_set_row_name (lp , ni+no+1 , "igual_a_uno");
    lpx_set_row_bnds (lp , ni+no+1 , LPX_FX , 1.0 , 1000.0 );

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no+1 );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    lpx_set_col_name (lp , nD+ni+no+1 , "te" );
    lpx_set_col_bnds (lp , nD+ni+no+1 , LPX_LO , 0.0 , 1000.0 );

    // completo la matriz de restricciones

    nf = ne;

    // primero meto la ultima columna

    for (s = 1 ; s <= ni ; s++ ){
      restri = *(xx + (i-1) + nD*(s-1) );
      if ( restri != 0.0 ){
	nf = nf + 1;
	ia[nf] = s;
	ja[nf] = nD+ni+no+1;
	ar[nf] = -restri;
      }
    }

    for (s = ni +  1 ; s <= ni+no ; s++ ){
      restri = *(yy + (i-1) + nD*(s-ni-1) );
      if ( restri != 0.0 ){
	nf = nf + 1;
	ia[nf] = s;
	ja[nf] = nD+ni+no+1;
	ar[nf] = -restri;
      }
    }

    nf = nf +1;
    ia[nf] = ni+no+1;
    ja[nf] = nD+ni+no+1;
    ar[nf] = 1;

    // ahora la ultima fila (salvo su ultimo elemento, que ya esta)

    for (s = 1 ; s <= no ; s++ ){
      restri = *(yy + (i-1) + nD*(s-1) );
      if ( restri != 0.0 ){
	nf = nf + 1;
	ia[nf] = ni+no+1;
	ja[nf] = nD+ni+s;
	ar[nf] = ( 1.0 / (no*restri) );
      }
    }


    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni ; j++){

      obje = *(xx + (i-1) + nD*(j-nD-1));

      if (obje != 0.0 ){
	lpx_set_obj_coef (lp , j , (-1.0/(ni*obje)));
      }
      //probablemente se puede poner cero el coeficiente si obje=0
    }


    for (j = nD+ni+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    lpx_set_obj_coef (lp , nD+ni+no+1 ,1.0 );

    // cargo la matriz de restricciones

    lpx_load_matrix (lp , nf , ia , ja , ar );

    // elimino mensajes coñazo
    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex
    rr = lpx_simplex ( lp );
    
    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion
    // el valor de la funcion objetivo ("ro") se va a la ultima columna

    *(zz + (i-1) + nD*(nD+ni+no) ) = lpx_get_obj_val(lp);

    // saco la "te"

    te = lpx_get_col_prim(lp , nD+ni+no+1);

    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) / te ;
    }

    // luego los slacks

    for(j = nD + 1 ; j <= nD+ni+no ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) / te ;
    }


    // borro el problema


    lpx_delete_prob(lp);
    


  }
  // fin del bucle gordo

} 

/*********************************************************/
/**************    fin de la funcion sbm      ************/
/*********************************************************/




/********************************************************/
/*                                                      */
/*       comienza el modelo slack based-measure         */
/*                     tipo BCC                         */
/*                                                      */
/********************************************************/


void sbm_bcc (char  ** nombre_DMUs,
	      int    * numero_DMUs,
	      char  ** nombre_inputs,
	      int    * numero_inputs,
	      char  ** nombre_outputs,
	      int    * numero_outputs,
	      double * xx,
	      double * yy,
	      int    * presolver,
	      int    * dual_simplex,
	      int    * information,
	      double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double te;
  double limi;

  int ne_max =( (ni + no + 1) * (nD+3));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no+2 filas y nD+ni+no+1 columnas
  
  int ne = 0;
  int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    /* añado la restricción de suma de lambdas igual a uno*/
    ne = ne + 1;
    ia[ne] = ni + no + 1;
    ja[ne] = j;
    ar[ne] = 1.0;

  }


    for (j = 1 ; j <= ni ; j++ ) {
      ne = ne+1;
      ia[ne] = j ;
      ja[ne] = j + nD ;
      ar[ne] = 1.0 ;
    }

    for (j = ni +  1 ; j <= ni+no ; j++ ) {
      ne = ne+1;
      ia[ne] = j ;
      ja[ne] = j + nD ;
      ar[ne] = -1.0 ;
    }


  
  // esta parte de la matriz de restricciones (la que corresponde a las
  // primeras nD+ni+no columnas) ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){


    // creo el problema, le asigno nombre e indico que es
    // de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_sbm_bcc");
    lpx_set_obj_dir (lp , LPX_MIN);

    // añado las filas, que van a ser ni+no+2, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no+2 );

    for ( t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , 0.0 , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , 0.0 , 1000.0 );
    }

    lpx_set_row_name (lp , ni+no+1 , "igual_a_uno_lambdas");
    lpx_set_row_bnds (lp , ni+no+1 , LPX_FX , 1.0 , 1000.0 );

    lpx_set_row_name (lp , ni+no+2 , "igual_a_uno");
    lpx_set_row_bnds (lp , ni+no+2 , LPX_FX , 1.0 , 1000.0 );

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no+1 );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    lpx_set_col_name (lp , nD+ni+no+1 , "te" );
    lpx_set_col_bnds (lp , nD+ni+no+1 , LPX_LO , 0.0 , 1000.0 );

    // completo la matriz de restricciones

    nf = ne;

    // primero meto la ultima columna, la que llamo "te"

    for (s = 1 ; s <= ni ; s++ ){
      restri = *(xx + (i-1) + nD*(s-1) );
      if ( restri != 0.0 ){
	nf = nf + 1;
	ia[nf] = s;
	ja[nf] = nD+ni+no+1;
	ar[nf] = -restri;
      }
    }

    for (s = ni +  1 ; s <= ni+no ; s++ ){
      restri = *(yy + (i-1) + nD*(s-ni-1) );
      if ( restri != 0.0 ){
	nf = nf + 1;
	ia[nf] = s;
	ja[nf] = nD+ni+no+1;
	ar[nf] = -restri;
      }
    }

    nf = nf +1;
    ia[nf] = ni+no+2;
    ja[nf] = nD+ni+no+1;
    ar[nf] = 1;

    // ahora la ultima fila (salvo su ultimo elemento, que ya esta)

    for (s = 1 ; s <= no ; s++ ){
      restri = *(yy + (i-1) + nD*(s-1) );
      if ( restri != 0.0 ){
	nf = nf + 1;
	ia[nf] = ni+no+2;
	ja[nf] = nD+ni+s;
	ar[nf] = ( 1.0 / (no*restri) );
      }
    }


    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni ; j++){

      obje = *(xx + (i-1) + nD*(j-nD-1));

      if (obje != 0.0 ){
	lpx_set_obj_coef (lp , j , (-1.0/(ni*obje)));
      }
      //probablemente se puede poner cero el coeficiente si obje=0
    }


    for (j = nD+ni+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    lpx_set_obj_coef (lp , nD+ni+no+1 ,1.0 );

    // cargo la matriz de restricciones
    lpx_load_matrix (lp , nf , ia , ja , ar );

    // elimino mensajes coñazo
    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex
    rr = lpx_simplex ( lp );
    

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion
    // el valor de la funcion objetivo ("ro") se va a la ultima columna

    *(zz + (i-1) + nD*(nD+ni+no) ) = lpx_get_obj_val(lp);

    // saco la "te"

    te = lpx_get_col_prim(lp , nD+ni+no+1);

    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) / te ;
    }

    // luego los slacks

    for(j = nD + 1 ; j <= nD+ni+no ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) / te ;
    }


    // borro el problema


    lpx_delete_prob(lp);
    
    /***********  acaba la segunda etapa   ***************/


  }
  // fin del bucle gordo

} 

/*********************************************************/
/************    fin de la funcion sbm_bcc    ************/
/*********************************************************/




/************************************************************/
/*                                                          */
/* comienza el modelo slack based-measure input-orientado   */
/*                                                          */
/*                             tipo CCR                     */
/************************************************************/


void sbm_ccr_io (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double te;
  double limi;

  int ne_max =( (ni + no) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no filas y nD+ni+no columnas
  
  int ne = 0;
  //int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }
  }


  for (j = 1 ; j <= ni ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = 1.0 ;
  }

  for (j = ni +  1 ; j <= ni+no ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = -1.0 ;
  }


  
  // esta parte de la matriz de restricciones
  // ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){


    // creo el problema, le asigno nombre e indico que es
    // de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_sbm_io");
    lpx_set_obj_dir (lp , LPX_MIN);

    // añado las filas, que van a ser ni+no, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no );

    for ( t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , *(xx + (i-1) + nD*(t-1)) , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , *(yy +(i-1)+ nD*(t-ni-1)) , 1000.0 );
    }


    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni ; j++){

      obje = *(xx + (i-1) + nD*(j-nD-1));

      if (obje != 0.0 ){
	lpx_set_obj_coef (lp , j , (-1.0/(ni*obje)));
      }
      //probablemente se puede poner cero el coeficiente si obje=0
    }


    for (j = nD+ni+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    lpx_set_obj_coef (lp , 0 , 1.0 );

    // cargo la matriz de restricciones
    lpx_load_matrix (lp , ne , ia , ja , ar );

    // elimino mensajes coñazo
    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex
    rr = lpx_simplex ( lp );
    

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion
    // el valor de la funcion objetivo ("ro.in") se va a la ultima columna

    *(zz + (i-1) + nD*(nD+ni+no) ) = lpx_get_obj_val(lp);


    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }

    // luego los slacks

    for(j = nD + 1 ; j <= nD+ni+no ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }


    // borro el problema


    lpx_delete_prob(lp);
    

  }
  // fin del bucle gordo

} 

/*********************************************************/
/************   fin de la funcion sbm_ccr_io  ************/
/*********************************************************/



/************************************************************/
/*                                                          */
/* comienza el modelo slack based-measure input-orientado   */
/*                                                          */
/*                  TIPO BCC                                */
/************************************************************/


  void sbm_bcc_io (char  ** nombre_DMUs,
                   int    * numero_DMUs,
 	           char  ** nombre_inputs,
	           int    * numero_inputs,
	           char  ** nombre_outputs,
	           int    * numero_outputs,
	           double * xx,
	           double * yy,
		   int    * presolver,
		   int    * dual_simplex,
		   int    * information,
	           double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double te;
  double limi;

  int ne_max =( (ni + no + 1) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no+1 filas y nD+ni+no columnas
  
  int ne = 0;
  //int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala; 

  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }
    
    /* añado la restricción  de suma de lambdas igual a uno */
    ne = ne +1;
    ia[ne] = ni + no +1;
    ja[ne] = j;
    ar[ne] = 1.0;
  }


  for (j = 1 ; j <= ni ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = 1.0 ;
  }

  for (j = ni +  1 ; j <= ni+no ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = -1.0 ;
  }


  
  // esta parte de la matriz de restricciones
  // ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){


    // creo el problema, le asigno nombre e indico que es
    // de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_sbm_io");
    lpx_set_obj_dir (lp , LPX_MIN);

    // añado las filas, que van a ser ni+no+1, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no+1 );

    for ( t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , *(xx + (i-1) + nD*(t-1)) , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , *(yy +(i-1)+ nD*(t-ni-1)) , 1000.0 );
    }

    lpx_set_row_name ( lp , ni+no+1 , "igual_a_uno" );
    lpx_set_row_bnds ( lp , ni+no+1 , LPX_FX , 1.0 , 1000.0 );

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni ; j++){

      obje = *(xx + (i-1) + nD*(j-nD-1));

      if (obje != 0.0 ){
	lpx_set_obj_coef (lp , j , (-1.0/(ni*obje)));
      }
      //probablemente se puede poner cero el coeficiente si obje=0
    }


    for (j = nD+ni+1 ; j <= nD+ni+no ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    lpx_set_obj_coef (lp , 0 , 1.0 );

    // cargo la matriz de restricciones
    lpx_load_matrix (lp , ne , ia , ja , ar );

    // elimino mensajes coñazo
    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex
    rr = lpx_simplex ( lp );
    

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion
    // el valor de la funcion objetivo ("ro.in") se va a la ultima columna

    *(zz + (i-1) + nD*(nD+ni+no) ) = lpx_get_obj_val(lp);


    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }

    // luego los slacks

    for(j = nD + 1 ; j <= nD+ni+no ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }


    // borro el problema


    lpx_delete_prob(lp);
    

  }
  // fin del bucle gordo

} 

/*********************************************************/
/************    fin de la funcion sbm_bcc_io   **********/
/*********************************************************/





/************************************************************/
/*                                                          */
/* comienza el modelo slack based-measure output-orientado  */
/*                                                          */
/*                      tipo CCR                            */
/************************************************************/


void sbm_ccr_oo (char  ** nombre_DMUs,
		 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double te;
  double limi;

  int ne_max =( (ni + no) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no filas y nD+ni+no columnas
  
  int ne = 0;
  //int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }
  }


  for (j = 1 ; j <= ni ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = 1.0 ;
  }

  for (j = ni +  1 ; j <= ni+no ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = -1.0 ;
  }


  
  // esta parte de la matriz de restricciones
  // ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){


    // creo el problema, le asigno nombre e indico que es
    // de maximización

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_sbm_oo");
    lpx_set_obj_dir (lp , LPX_MAX);

    // añado las filas, que van a ser ni+no, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no );

    for ( t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , *(xx + (i-1) + nD*(t-1)) , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , *(yy +(i-1)+ nD*(t-ni-1)) , 1000.0 );
    }


    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni ; j++){
      lpx_set_obj_coef (lp , j , 0.0); 
    }


    for (j = nD+ni+1 ; j <= nD+ni+no ; j++ ) {
      obje = *(yy + (i-1) + nD*(j-nD-ni-1));
      if (obje != 0.0) {
        lpx_set_obj_coef (lp , j , ( 1.0/( no*obje)));
      }
    }

    lpx_set_obj_coef (lp , 0 , 1.0 ); 

    // cargo la matriz de restricciones
    lpx_load_matrix (lp , ne , ia , ja , ar );

    // elimino mensajes coñazo
    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex
    rr = lpx_simplex ( lp );
    

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 
    
    // extraigo la solucion
    // el valor de la funcion objetivo ("ro.out") se va a la ultima columna
    // luego tendré que invertirlo

    *(zz + (i-1) + nD*(nD+ni+no) ) = lpx_get_obj_val(lp);


    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }

    // luego los slacks

    for(j = nD + 1 ; j <= nD+ni+no ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }


    // borro el problema


    lpx_delete_prob(lp);
    

  }
  // fin del bucle gordo

} 

/*********************************************************/
/************   fin de la funcion sbm_ccr_oo  ************/
/*********************************************************/



/************************************************************/
/*                                                          */
/* comienza el modelo slack based-measure output-orientado  */
/*                                                          */
/*                  TIPO BCC                                */
/************************************************************/


void sbm_bcc_oo (char  ** nombre_DMUs,
                 int    * numero_DMUs,
		 char  ** nombre_inputs,
		 int    * numero_inputs,
		 char  ** nombre_outputs,
		 int    * numero_outputs,
		 double * xx,
		 double * yy,
		 int    * presolver,
		 int    * dual_simplex,
		 int    * information,
		 double * zz) {


  int nD = * numero_DMUs;
  int ni = * numero_inputs;
  int no = * numero_outputs;

  int PSLV = * presolver;
  int DUAL = * dual_simplex;
  int INFR = * information;

/*
  Rprintf("numero de DMUs: %d \n",nD);
  Rprintf("numero de inputs: %d \n",ni);
  Rprintf("numero de outputs: %d \n",no);
*/

  int i, j, k, s,t;

  char  tt[] = "L_";
  char  ss[256] ;

  strcpy(ss, tt);

  char smas[256];
  char smenos[256];
  char tmas[] = "s+";
  char tmenos[] = "s-";

  strcpy (smas, tmas);
  strcpy (smenos , tmenos);

  // obje lo uso para rellenar los coeficientes de la funcion objetivo
  // y restri para los coeficientes de la matriz de restricciones.

  double obje, restri;
  double te;
  double limi;

  int ne_max =( (ni + no + 1) * (nD+1));
  // contiene el numero de elementos de la matriz
  // de restricciones; algunos pueden ser nulos.
  // esta matriz tiene ni+no+1 filas y nD+ni+no columnas
  
  int ne = 0;
  //int nf;

  int ia[1 + ne_max];
  int ja[1 + ne_max];
  double ar[1 + ne_max];

  LPX * lp;

  LPXKKT  kkt;
  int rr, scala;
 
  // establecemos la matriz de restricciones rellenando los 
  // arrays ia, ja, ar


  for ( j = 1 ; j <= nD ; j++ ) {

    for ( i = 1 ; i <= ni ; i++ ) {

      restri = *( xx + (j-1) +  nD*(i-1));

      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }

    for ( i = ni+1 ; i <= ni+no ; i++){

      restri = *( yy + (j-1) + nD*(i-ni-1));
      
      if (restri != 0.0 ) {
	ne = ne+1;
	ia[ne] = i;
	ja[ne] = j;
	ar[ne] = restri;
      }
    }
    
    /* añado la restricción  de suma de lambdas igual a uno */
    ne = ne +1;
    ia[ne] = ni + no +1;
    ja[ne] = j;
    ar[ne] = 1.0;
  }


  for (j = 1 ; j <= ni ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = 1.0 ;
  }

  for (j = ni +  1 ; j <= ni+no ; j++ ) {
    ne = ne+1;
    ia[ne] = j ;
    ja[ne] = j + nD ;
    ar[ne] = -1.0 ;
  }


  
  // esta parte de la matriz de restricciones
  // ya no va a cambiar. El numero de elementos no nulos
  // que tiene es ne.


    /*********     BUCLE GORDO      *********/


  // ahora hago un bucle gordo sobre las DMUs
  // porque para cada una de ellas habra un problema


  for ( i = 1 ; i <= nD ; i++ ){


    // creo el problema, le asigno nombre e indico que es
    // de minimizacion

    lp = lpx_create_prob();
    lpx_set_prob_name (lp , "dea_sbm_bcc_oo");
    lpx_set_obj_dir (lp , LPX_MAX);

    // añado las filas, que van a ser ni+no+1, les doy nombre
    // y establezco sus limites

    lpx_add_rows (lp ,  ni+no+1 );

    for ( t = 1 ; t <= ni ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_inputs + (t-1)));
      lpx_set_row_bnds (lp , t , LPX_FX , *(xx + (i-1) + nD*(t-1)) , 1000.0 );
    }

    for (t = ni+1 ; t <= ni+no ; t++ ) {

      lpx_set_row_name (lp , t , *(nombre_outputs + (t-ni-1) ) );
      lpx_set_row_bnds (lp , t , LPX_FX , *(yy +(i-1)+ nD*(t-ni-1)) , 1000.0 );
    }

    lpx_set_row_name ( lp , ni+no+1 , "igual_a_uno" );
    lpx_set_row_bnds ( lp , ni+no+1 , LPX_FX , 1.0 , 1000.0 );

    // ahora añado columnas, es decir variables estructurales
    // añado tambien sus nombres y sus limites

    lpx_add_cols (lp , nD+ni+no );

    for (t = 1 ; t <= nD ; t++){
      lpx_set_col_name (lp , t , strcat(ss, *(nombre_DMUs + (t-1) ) ) );
      strcpy(ss , tt );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+1 ; t <= nD+ni ; t++ ){
      lpx_set_col_name(lp , t , strcat (smenos , *(nombre_inputs + (t-nD-1))));
      strcpy(smenos , tmenos );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }

    for (t = nD+ni+1 ; t <= nD+ni+no ; t++ ) {
      lpx_set_col_name(lp ,t , strcat(smas , *(nombre_outputs + (t-nD-ni-1))));
      strcpy (smas , tmas );
      lpx_set_col_bnds (lp , t , LPX_LO , 0.0 , 1000.0 );
    }


    // establezco la funcion objetivo

    for ( j = 1 ; j <= nD ; j++ ) {
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+1 ; j <= nD+ni ; j++){
      lpx_set_obj_coef (lp , j , 0.0 );
    }

    for (j = nD+ni+1 ; j <= nD+ni+no ; j++ ) {
      obje = *(yy + (i-1) + nD*(j-nD-ni-1));

      if (obje != 0.0 ){
	lpx_set_obj_coef (lp , j , (1.0/(no*obje)));
      }
      //probablemente se puede poner cero el coeficiente si obje=0
     }
    

    lpx_set_obj_coef (lp , 0 , 1.0 );

    // cargo la matriz de restricciones
    lpx_load_matrix (lp , ne , ia , ja , ar );

    // elimino mensajes coñazo
    lpx_set_int_parm (lp , LPX_K_MSGLEV , 1 );

    // flag del presolver
    lpx_set_int_parm ( lp , LPX_K_PRESOL , PSLV);

    // dual simplex
    lpx_set_int_parm ( lp , LPX_K_DUAL , DUAL );

    // aplico el simplex
    rr = lpx_simplex ( lp );
    

    if ( INFR != 0 ) {
      if ( lpx_get_status ( lp ) == LPX_OPT ) { 
	lpx_check_kkt ( lp , 1 , &kkt );
	Rprintf ("problem %d: OK, %c %c %c %c \n", i , kkt.pe_quality , kkt.pb_quality , kkt.de_quality , kkt.db_quality );
      }
      else
	Rprintf ("problem %d: FAULT \n", i );
    }; 

    // extraigo la solucion
    // el valor de la funcion objetivo ("ro.out") se va a la ultima columna

    *(zz + (i-1) + nD*(nD+ni+no) ) = lpx_get_obj_val(lp);


    // el valor de las lambdas se va a las primeras nD columnas

    for(j = 1 ; j <= nD ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }

    // luego los slacks

    for(j = nD + 1 ; j <= nD+ni+no ; j++) {
      *(zz + (i-1) + (j-1)*nD ) = lpx_get_col_prim(lp , j) ;
    }


    // borro el problema


    lpx_delete_prob(lp);
    

  }
  // fin del bucle gordo

} 

/*********************************************************/
/************    fin de la funcion sbm_bcc_oo   **********/
/*********************************************************/
