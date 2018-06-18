/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Chuanfu Luo (luochuanfu@gmail.com)
   Contributing author: Luciano Laratelli (Luciano@Laratelli.com)
   Contributing author: Preston Moore (P.Moore@Usciences.edu)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "angle_table3D.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum{LINEAR,SPLINE};

#define MAXLINE 1024
#define SMALL 0.001
#define TINY  1.E-10

/* ---------------------------------------------------------------------- */

AngleTable3D::AngleTable3D(LAMMPS *lmp) : Angle(lmp)
{
    writedata = 0;
    ntables = 0;
    tables = nullptr;
}

/* ---------------------------------------------------------------------- */

AngleTable3D::~AngleTable3D()
{
    for (int m = 0; m < ntables; m++) free_table(&tables[m]);
    memory->sfree(tables);

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(theta0);
        memory->destroy(tabindex);
    }
}

/* ---------------------------------------------------------------------- */

void AngleTable3D::compute(int eflag, int vflag)
{
    int i1,i2,i3,n,type;
    double eangle,f1[3],f3[3];
    double delx1,dely1,delz1,delx2,dely2,delz2;
    double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;
    double theta,u,mdu_theta, mdu_r1, mdu_r2; //mdu: minus du, -du/dx=f

    eangle = 0.0;
    if (eflag || vflag) ev_setup(eflag,vflag);
    else evflag = 0;

    double **x = atom->x;
    double **f = atom->f;
    int **anglelist = neighbor->anglelist;
    int nanglelist = neighbor->nanglelist;
    int nlocal = atom->nlocal;
    int newton_bond = force->newton_bond;

    for (n = 0; n < nanglelist; n++) {
        i1 = anglelist[n][0];
        i2 = anglelist[n][1];
        i3 = anglelist[n][2];
        type = anglelist[n][3];

        // 1st bond

        delx1 = x[i1][0] - x[i2][0];
        dely1 = x[i1][1] - x[i2][1];
        delz1 = x[i1][2] - x[i2][2];

        rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
        r1 = sqrt(rsq1);

        // 2nd bond

        delx2 = x[i3][0] - x[i2][0];
        dely2 = x[i3][1] - x[i2][1];
        delz2 = x[i3][2] - x[i2][2];

        rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
        r2 = sqrt(rsq2);

        // angle (cos and sin)

        c = delx1*delx2 + dely1*dely2 + delz1*delz2;
        c /= r1*r2;

        if (c > 1.0) c = 1.0;
        if (c < -1.0) c = -1.0;

        s = sqrt(1.0 - c*c);
        if (s < SMALL) s = SMALL;
        s = 1.0/s;

        // tabulated force & energy

        theta = acos(c);
        uf_lookup(type,theta,r1, r2, u, mdu_theta, mdu_r1, mdu_r2);

        if (eflag) eangle = u;


        a = mdu_theta * s;
        a11 = a*c / rsq1;
        a12 = -a / (r1*r2);
        a22 = a*c / rsq2;

        f1[0] = a11*delx1 + a12*delx2;
        f1[1] = a11*dely1 + a12*dely2;
        f1[2] = a11*delz1 + a12*delz2;
        f3[0] = a22*delx2 + a12*delx1;
        f3[1] = a22*dely2 + a12*dely1;
        f3[2] = a22*delz2 + a12*delz1;

        // apply force to each of 3 atoms

        if (newton_bond || i1 < nlocal) {
            f[i1][0] += f1[0];
            f[i1][1] += f1[1];
            f[i1][2] += f1[2];
        }

        if (newton_bond || i2 < nlocal) {
            f[i2][0] -= f1[0] + f3[0];
            f[i2][1] -= f1[1] + f3[1];
            f[i2][2] -= f1[2] + f3[2];
        }

        if (newton_bond || i3 < nlocal) {
            f[i3][0] += f3[0];
            f[i3][1] += f3[1];
            f[i3][2] += f3[2];
        }

        double fbond_r1 = mdu_r1/r1;

        // apply force to each of 2 atoms

        if (newton_bond || i1 < nlocal) {
            f[i1][0] -= delx1* fbond_r1;
            f[i1][1] -= dely1* fbond_r1;
            f[i1][2] -= delz1* fbond_r1;
        }

        if (newton_bond || i2 < nlocal) {
            f[i2][0] += delx1*fbond_r1;
            f[i2][1] += dely1*fbond_r1;
            f[i2][2] += delz1*fbond_r1;
        }

        double fbond_r2 = mdu_r2/r2;

        // apply force to each of 2 atoms

        if (newton_bond || i1 < nlocal) {
            f[i2][0] += delx2* fbond_r2;
            f[i2][1] += dely2* fbond_r2;
            f[i2][2] += delz2* fbond_r2;
        }

        if (newton_bond || i2 < nlocal) {
            f[i3][0] -= delx2*fbond_r2;
            f[i3][1] -= dely2*fbond_r2;
            f[i3][2] -= delz2*fbond_r2;
        }

        if (evflag) ev_tally(i1,i2,i3,nlocal,newton_bond,eangle,f1,f3,
                             delx1,dely1,delz1,delx2,dely2,delz2);
    }
}

/* ---------------------------------------------------------------------- */

void AngleTable3D::allocate()
{
    allocated = 1;
    int n = atom->nangletypes;

    memory->create(theta0,n+1,"angle:theta0");
    memory->create(tabindex,n+1,"angle:tabindex");

    memory->create(setflag,n+1,"angle:setflag");
    for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void AngleTable3D::settings(int narg, char **arg)
{
    if (narg != 2) error->all(FLERR,"Illegal angle_style command");

    if (strcmp(arg[0],"linear") == 0) tabstyle = LINEAR;
    else if (strcmp(arg[0],"spline") == 0) tabstyle = SPLINE;
    else error->all(FLERR,"Unknown table style in angle style table");

    tablength = force->inumeric(FLERR,arg[1]);
    if (tablength < 2) error->all(FLERR,"Illegal number of angle table entries");

    // delete old tables, since cannot just change settings

    for (int m = 0; m < ntables; m++) free_table(&tables[m]);
    memory->sfree(tables);

    if (allocated) {
        memory->destroy(setflag);
        memory->destroy(tabindex);
    }
    allocated = 0;

    ntables = 0;
    tables = NULL;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

//AKA MAIN LOOP
void AngleTable3D::coeff(int narg, char **arg)
{
    if (narg != 3) error->all(FLERR,"Illegal angle_coeff command");
    if (!allocated) allocate();

    int ilo,ihi;
    force->bounds(FLERR,arg[0],atom->nangletypes,ilo,ihi);

    int me;
    MPI_Comm_rank(world,&me);
    tables = (Table *)
        memory->srealloc(tables,(ntables+1)*sizeof(Table),"angle:tables");
    Table *tb = &tables[ntables];
    null_table(tb);
    if (me == 0) read_table(tb,arg[1],arg[2]);
    bcast_table(tb);

    // error check on table parameters
    if(tb->angle_low < 0.0 || tb->angle_high > 180.0)
    {
        error->one(FLERR, "Angle must be between 0.0 and 180.0");
    } else if (tb->r1_low < TINY || tb->r2_low < TINY) {
        error->one(FLERR, "Minimum r distance is too small");
    }
    //end error check
    if (tb->input_count <= 1) error->one(FLERR,"Invalid combined table length (ntheta * nr1 * nr2 * 4)");

    spline_table(tb);
    compute_table(tb);

    // store ptr to table in tabindex

    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
        tabindex[i] = ntables;
        setflag[i] = 1;
        theta0[i] = tb->theta0;
        count++;
    }
    ntables++;

    if (count == 0) error->all(FLERR,"Illegal angle_coeff command");
}

/* ----------------------------------------------------------------------
   return an equilbrium angle length
   should not be used, since don't know minimum of tabulated function
------------------------------------------------------------------------- */

double AngleTable3D::equilibrium_angle(int i)
{
    return theta0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
 ------------------------------------------------------------------------- */

void AngleTable3D::write_restart(FILE *fp)
{
    fwrite(&tabstyle,sizeof(int),1,fp);
    fwrite(&tablength,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
    proc 0 reads from restart file, bcasts
 ------------------------------------------------------------------------- */

void AngleTable3D::read_restart(FILE *fp)
{
    if (comm->me == 0) {
        fread(&tabstyle,sizeof(int),1,fp);
        fread(&tablength,sizeof(int),1,fp);
    }
    MPI_Bcast(&tabstyle,1,MPI_INT,0,world);
    MPI_Bcast(&tablength,1,MPI_INT,0,world);

    allocate();
}

/* ---------------------------------------------------------------------- */

double AngleTable3D::single(int type, int i1, int i2, int i3)
{
    double **x = atom->x;

    double delx1 = x[i1][0] - x[i2][0];
    double dely1 = x[i1][1] - x[i2][1];
    double delz1 = x[i1][2] - x[i2][2];
    domain->minimum_image(delx1,dely1,delz1);
    double r1 = sqrt(delx1*delx1 + dely1*dely1 + delz1*delz1);

    double delx2 = x[i3][0] - x[i2][0];
    double dely2 = x[i3][1] - x[i2][1];
    double delz2 = x[i3][2] - x[i2][2];
    domain->minimum_image(delx2,dely2,delz2);
    double r2 = sqrt(delx2*delx2 + dely2*dely2 + delz2*delz2);

    double c = delx1*delx2 + dely1*dely2 + delz1*delz2;
    c /= r1*r2;
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    double theta = acos(c);
    double u = 0.0;
    double dummy = 0.0;
    double dummy2 = 0.0;
    double dummy3 = 0.0;
    uf_lookup(type,theta, r1, r2, u, dummy, dummy2, dummy3);
    return u;
}

/* ---------------------------------------------------------------------- */

void AngleTable3D::null_table(Table *tb)
{
    tb->e2file = tb->f2file = NULL;
    tb->ang = tb->e = tb->de = NULL;
}

/* ---------------------------------------------------------------------- */

void AngleTable3D::free_table(Table *tb)
{
    memory->destroy(tb->e2file);
    memory->destroy(tb->f2file);

    memory->destroy(tb->ang);
    memory->destroy(tb->e);
    memory->destroy(tb->de);
    memory->destroy(tb->df);
    memory->destroy(tb->e2);
    memory->destroy(tb->f2);
    free(tb->energies_forces);
}

/* ----------------------------------------------------------------------
   read table file, only called by proc 0
------------------------------------------------------------------------- */

void AngleTable3D::read_table(Table *tb, char *file, char *keyword)
{
    char line[MAXLINE];

    // open file

    FILE *fp = force->open_potential(file);
    if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open file %s",file);
        error->one(FLERR,str);
    }


    // loop until section found with matching keyword

    while (1) {
        if (fgets(line,MAXLINE,fp) == NULL)
            error->one(FLERR,"Did not find keyword in table file");
        if (strspn(line," \t\n") == strlen(line)) continue;    // blank line
        if (line[0] == '#') continue;                          // comment
        char *word = strtok(line," \t\n\r");
        if (strcmp(word,keyword) == 0) break;           // matching keyword
        else {
            error->one(FLERR, "YOU'RE A BIG DUMMY || THIS CODE SHOULD NOT EXECUTE");
        }
        /*
          fgets(line,MAXLINE,fp);                         // no match, skip section
          param_extract(tb,line);
          fgets(line,MAXLINE,fp);
          for (int i = 0; i < (tb->input_count/4); i++) fgets(line,MAXLINE,fp);
         */
    }
    fgets(line,MAXLINE,fp);
    param_extract(tb,line);

    tb->input_count = tb->ninput_theta * tb->ninput_r2 * tb->ninput_r1 * 4;

    tb->energies_forces = (double*)(malloc(sizeof(double*) * tb->input_count)) ;

    int theta_offset = tb->ninput_r1*tb->ninput_r2*4;
    int r1_offset = tb->ninput_r2*4;
    int r2_offset = 4;

    fgets(line,MAXLINE,fp);
    param_extract(tb,line);

    double temp_a, temp_r1, temp_r2;

    int itmp;
    for(int i = 0; i < tb->ninput_theta; i++)
    {
        for(int j = 0; j < tb->ninput_r1; j++)
        {
            for(int k = 0; k < tb->ninput_r2; k++)
            {
                fgets(line, MAXLINE, fp);
                int index = i*theta_offset + j*r1_offset + k*r2_offset;
                double &u = tb->energies_forces[index + 0];
                double &f_theta = tb->energies_forces[index + 1];
                double &f_r1 = tb->energies_forces[index + 2];
                double &f_r2 = tb->energies_forces[index + 3];

                sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg",
                       &itmp, &temp_a, &temp_r1, &temp_r2, &u, &f_theta, &f_r1, &f_r2);
                //assert(index == (itmp*4));
                assert(index == itmp);
                if(index == 0)
                {
                    assert(temp_a == tb->angle_low);
                    assert(temp_r1 == tb->r1_low);
                    assert(temp_r2 == tb->r2_low);
                }
                else if(index == (tb->ninput_theta * tb->ninput_r1 * tb->ninput_r2 * 4))
                {
                    assert(temp_a == tb->angle_high);
                    assert(temp_r1 == tb->r1_high);
                    assert(temp_r2 == tb->r2_high);
                }
                //std::cout << index << std::endl;
                //std::cout << line << std::endl;
                //std::cout << i << " " << j << " " << k << std::endl;
            }
        }
    }

}

/* ----------------------------------------------------------------------
   build spline representation of e,f over entire range of read-in table
   this function sets these values in e2file,f2file
------------------------------------------------------------------------- */

void AngleTable3D::spline_table(Table *tb)
{
    /*
  memory->create(tb->e2file,tb->ninput,"angle:e2file");
  memory->create(tb->f2file,tb->ninput,"angle:f2file");

  double ep0 = - tb->fafile[0];
  double epn = - tb->fafile[tb->ninput-1];
  spline(tb->afile,tb->efile,tb->ninput,ep0,epn,tb->e2file);

  if (tb->fpflag == 0) {
    tb->fplo = (tb->fafile[1] - tb->fafile[0]) / (tb->afile[1] - tb->afile[0]);
    tb->fphi = (tb->fafile[tb->ninput-1] - tb->fafile[tb->ninput-2]) /
        (tb->afile[tb->ninput-1] - tb->afile[tb->ninput-2]);
  }

  double fp0 = tb->fplo;
  double fpn = tb->fphi;
  spline(tb->afile,tb->fafile,tb->ninput,fp0,fpn,tb->f2file);
     */
}

/* ----------------------------------------------------------------------
   compute a,e,f vectors from splined values
------------------------------------------------------------------------- */

void AngleTable3D::compute_table(Table *tb)
{
    // delta = table spacing in angle for N-1 bins

    tb->delta_theta = (tb->angle_high - tb->angle_low) / (tb->ninput_theta - 1);
    tb->delta_r1 = (tb->r1_high - tb->r1_low) / (tb->ninput_r1 - 1);
    tb->delta_r2 = (tb->r2_high - tb->r2_low) / (tb->ninput_r2 - 1);

    tb->inv_delta_theta = 1.0 / tb->delta_theta;
    tb->inv_delta_r1 = 1.0 / tb->delta_r1;
    tb->inv_delta_r2 = 1.0 / tb->delta_r2;

    // N-1 evenly spaced bins in angle from 0 to PI
    // ang,e,f = value at lower edge of bin
    // de,df values = delta values of e,f
    // ang,e,f are N in length so de,df arrays can compute difference
}

/* ----------------------------------------------------------------------
   extract attributes from parameter line in table section
   format of line: N value FP fplo fphi EQ theta0
   N is required, other params are optional
------------------------------------------------------------------------- */

void AngleTable3D::param_extract(Table *tb, char *line)
{
    tb->ninput = 0;
    tb->fpflag = 0;
    tb->theta0 = 180.0;

    char *word = strtok(line," \t\n\r\f");
    while (word) {
        if (strcmp(word,"THETA") == 0) {
            word = strtok(NULL," \t");
            tb->ninput_theta = atoi(word);
            word = strtok(NULL," \t");
            tb->angle_low = atof(word);
            word = strtok(NULL," \t");
            tb->angle_high = atof(word);
        } else if (strcmp(word,"R1") == 0) {
            word = strtok(NULL," \t");
            tb->ninput_r1 = atoi(word);
            word = strtok(NULL," \t");
            tb->r1_low = atof(word);
            word = strtok(NULL," \t");
            tb->r1_high = atof(word);
        } else if (strcmp(word,"R2") == 0) {
            word = strtok(NULL," \t");
            tb->ninput_r2 = atoi(word);
            word = strtok(NULL," \t");
            tb->r2_low = atof(word);
            word = strtok(NULL," \t");
            tb->r2_high = atof(word);
        } else {
            error->one(FLERR,"Invalid keyword in angle table 3D parameters");
        }
        word = strtok(NULL," \t\n\r\f");
    }
    //if (tb->ninput == 0) error->one(FLE RR,"Angle table parameters did not set N");
}

/* ----------------------------------------------------------------------
   broadcast read-in table info from proc 0 to other procs
   this function communicates these values in Table:
     ninput,afile,efile,fafile,fpflag,fplo,fphi,theta0
------------------------------------------------------------------------- */

void AngleTable3D::bcast_table(Table *tb)
{
    MPI_Bcast(&tb->ninput,1,MPI_INT,0,world);

    int me;
    MPI_Comm_rank(world,&me);
    if (me > 0) {
        memory->create(tb->energies_forces,tb->input_count,"angle:energies_forces");
    }

    MPI_Bcast(tb->energies_forces,tb->input_count,MPI_DOUBLE,0,world);

    MPI_Bcast(&tb->fpflag,1,MPI_INT,0,world);
    if (tb->fpflag) {
        MPI_Bcast(&tb->fplo,1,MPI_DOUBLE,0,world);
        MPI_Bcast(&tb->fphi,1,MPI_DOUBLE,0,world);
    }
    MPI_Bcast(&tb->theta0,1,MPI_DOUBLE,0,world);
}

/* ----------------------------------------------------------------------
   spline and splint routines modified from Numerical Recipes
------------------------------------------------------------------------- */

/* "this routine stores an array, y2, with the second derivatives
 * of the interpolated function at the interpolated points, using
 * functional values" --Preston Moore, PhD
 */
void AngleTable3D::spline(double *x, double *y, int n,
                          double yp1, double ypn, double *y2)
{
    int i,k;
    double p,qn,sig,un;
    double *u = new double[n];

    if (yp1 > 0.99e30) y2[0] = u[0] = 0.0;
    else {
        y2[0] = -0.5;
        u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
    }
    for (i = 1; i < n-1; i++) {
        sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
        p = sig*y2[i-1] + 2.0;
        y2[i] = (sig-1.0) / p;
        u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
        u[i] = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
    }
    if (ypn > 0.99e30) qn = un = 0.0;
    else {
        qn = 0.5;
        un = (3.0/(x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
    }
    y2[n-1] = (un-qn*u[n-2]) / (qn*y2[n-2] + 1.0);
    for (k = n-2; k >= 0; k--) y2[k] = y2[k]*y2[k+1] + u[k];

    delete [] u;
}

/* ---------------------------------------------------------------------- */

double AngleTable3D::splint(double *xa, double *ya, double *y2a, int n, double x)
{
    int klo,khi,k;
    double h,b,a,y;

    klo = 0;
    khi = n-1;
    while (khi-klo > 1) {
        k = (khi+klo) >> 1;//why god
        if (xa[k] > x) khi = k;
        else klo = k;
    }
    h = xa[khi]-xa[klo];
    a = (xa[khi]-x) / h;
    b = (x-xa[klo]) / h;
    y = a*ya[klo] + b*ya[khi] +
        ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * (h*h)/6.0;
    return y;
}

/* ----------------------------------------------------------------------
   calculate potential u and force f at angle x
------------------------------------------------------------------------- */

void AngleTable3D::uf_lookup(int type, double theta, double r1, double r2, double &u, double &mdu_theta, double &mdu_r1, double &mdu_r2)
{
    if (!ISFINITE(theta)) {
        error->one(FLERR,"Illegal angle in angle style table");
    }

    double fraction_theta, fraction_r1, fraction_r2;
    const Table *tb = &tables[tabindex[type]];
    int itable_theta = static_cast<int> ((theta - tb->angle_low) * tb->inv_delta_theta);
    int itable_r1 = static_cast<int> ((r1 - tb->r1_low) * tb->inv_delta_r1);
    int itable_r2 = static_cast<int> ((r2 - tb->r2_low) * tb->inv_delta_r2 );

    /***** check boundaries for theta ******/
    if (itable_theta < 0)
    {
     itable_theta = 0;
     theta = tb->angle_low;
     std::cerr << "WARNING: Reached bottom of table in theta, rolling itable_theta up to 0" << std::endl;
    }
    else if (itable_theta >= tb->ninput_theta)
    {
        itable_theta = tb->ninput_theta - 1;
        theta = itable_theta * tb->delta_theta + tb->angle_low;
        std::cerr << "WARNING: Reached end of table in theta, rolling back to max_theta - 1" << std::endl;
    }
    /***** check boundaries for r1 ******/
    if (itable_r1 < 0)
    {
        itable_r1 = 0;
        r1 = tb->r1_low;
        std::cerr << "WARNING: Reached bottom of table in R1, rolling itable_r1 up to 0" << std::endl;
    }
    else if (itable_r1 >= tb->ninput_r1 )
    {
        itable_r1 = tb->ninput_r1 - 1;
        r1 = itable_r1 * tb->delta_r1 + tb->r1_low;
        std::cerr << "WARNING: Reached end of table in r1, rolling back to max_r1 - 1" << std::endl;
    }
    /***** check boundaries for r2 ******/
    if (itable_r2 < 0)
    {
        itable_r2 = 0;
        r2 = tb->r2_low;
        std::cerr << "WARNING: Reached bottom of table in R2, rolling itable_r2 up to 0" << std::endl;
    }
    else if (itable_r2 >= tb->ninput_r2 )
    {
        itable_r2 = tb->ninput_r2 - 1;
        r2 = itable_r2 * tb->delta_r2 + tb->r2_low;
        std::cerr << "WARNING: Reached end of table in r2, rolling back to max_r2-1" << std::endl;
    }

    using namespace std;

    if (tabstyle == LINEAR) {
        fraction_theta = (theta - (tb->angle_low + itable_theta * tb->delta_theta)) * tb->inv_delta_theta;
        fraction_r1 = (r1 - (tb->r1_low + itable_r1 * tb->delta_r1)) * tb->inv_delta_r1;
        fraction_r2 = (r2 - (tb->r2_low + itable_r2 * tb->delta_r2)) * tb->inv_delta_r2;

        /*
        cout << fraction_theta << " " << theta << " " << tb->angle_low + itable_theta * tb->delta_theta << endl;
        cout << fraction_r1 << " " << r1 << " " << tb->r1_low + itable_r1 * tb->delta_r1 << endl;
        cout << fraction_r2 << " " << r2 << " " << tb->r2_low + itable_r2 * tb->delta_r2 << endl;
        exit(10000);
         */

        int index_000 = (itable_theta)     * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1)     * tb-> ninput_r2 * 4 + (itable_r2)     * 4;
        int index_001 = (itable_theta)     * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1)     * tb-> ninput_r2 * 4 + (itable_r2 + 1) * 4;
        int index_010 = (itable_theta)     * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1 + 1) * tb-> ninput_r2 * 4 + (itable_r2)     * 4;
        int index_011 = (itable_theta)     * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1 + 1) * tb-> ninput_r2 * 4 + (itable_r2 + 1) * 4;
        int index_100 = (itable_theta + 1) * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1)     * tb-> ninput_r2 * 4 + (itable_r2)     * 4;
        int index_101 = (itable_theta + 1) * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1)     * tb-> ninput_r2 * 4 + (itable_r2 + 1) * 4;
        int index_110 = (itable_theta + 1) * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1 + 1) * tb-> ninput_r2 * 4 + (itable_r2)     * 4;
        int index_111 = (itable_theta + 1) * tb->ninput_r1 * tb->ninput_r2 * 4 + (itable_r1 + 1) * tb-> ninput_r2 * 4 + (itable_r2 + 1) * 4;

        u = 0.0;
        u+= (1-fraction_theta) * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_000];
        u+= (1-fraction_theta) * (1-fraction_r1) * (fraction_r2)   * tb->energies_forces[index_001];
        u+= (1-fraction_theta) * (fraction_r1)   * (1-fraction_r2) * tb->energies_forces[index_010];
        u+= (1-fraction_theta) * (fraction_r1)   * (fraction_r2)   * tb->energies_forces[index_011];
        u+= (fraction_theta)   * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_100];
        u+= (fraction_theta)   * (1-fraction_r1) * (fraction_r2)   * tb->energies_forces[index_101];
        u+= (fraction_theta)   * (fraction_r1)   * (1-fraction_r2) * tb->energies_forces[index_110];
        u+= (fraction_theta)   * (fraction_r1)   * (fraction_r2)   * tb->energies_forces[index_111];

        mdu_theta = 0.0;
        mdu_theta+= (1-fraction_theta) * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_000 + 1];
        mdu_theta+= (1-fraction_theta) * (1-fraction_r1) * (fraction_r2) * tb->energies_forces[index_001 + 1];
        mdu_theta+= (1-fraction_theta) * (fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_010 + 1];
        mdu_theta+= (1-fraction_theta) * (fraction_r1) * (fraction_r2) * tb->energies_forces[index_011 + 1];
        mdu_theta+= (fraction_theta) * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_100 + 1];
        mdu_theta+= (fraction_theta) * (1-fraction_r1) * (fraction_r2) * tb->energies_forces[index_101 + 1];
        mdu_theta+= (fraction_theta) * (fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_110 + 1];
        mdu_theta+= (fraction_theta) * (fraction_r1) * (fraction_r2) * tb->energies_forces[index_111 + 1];

        mdu_r1 = 0.0;
        mdu_r1+= (1-fraction_theta) * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_000 + 2];
        mdu_r1+= (1-fraction_theta) * (1-fraction_r1) * (fraction_r2) * tb->energies_forces[index_001 + 2];
        mdu_r1+= (1-fraction_theta) * (fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_010 +2];
        mdu_r1+= (1-fraction_theta) * (fraction_r1) * (fraction_r2) * tb->energies_forces[index_011 +2];
        mdu_r1+= (fraction_theta) * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_100 +2];
        mdu_r1+= (fraction_theta) * (1-fraction_r1) * (fraction_r2) * tb->energies_forces[index_101 +2];
        mdu_r1+= (fraction_theta) * (fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_110 +2];
        mdu_r1+= (fraction_theta) * (fraction_r1) * (fraction_r2) * tb->energies_forces[index_111 +2];

        mdu_r2 = 0.0;
        mdu_r2+= (1-fraction_theta) * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_000 +3];
        mdu_r2+= (1-fraction_theta) * (1-fraction_r1) * (fraction_r2) * tb->energies_forces[index_001 +3];
        mdu_r2+= (1-fraction_theta) * (fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_010 +3];
        mdu_r2+= (1-fraction_theta) * (fraction_r1) * (fraction_r2) * tb->energies_forces[index_011 +3];
        mdu_r2+= (fraction_theta) * (1-fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_100 +3];
        mdu_r2+= (fraction_theta) * (1-fraction_r1) * (fraction_r2) * tb->energies_forces[index_101 +3];
        mdu_r2+= (fraction_theta) * (fraction_r1) * (1-fraction_r2) * tb->energies_forces[index_110 +3];
        mdu_r2+= (fraction_theta) * (fraction_r1) * (fraction_r2) * tb->energies_forces[index_111 +3];

    } else {
        error->one(FLERR,"ILLEGAL METHOD OF INTERPOLATION: LINEAR-ONLY FOR NOW BUB");
    }
    /*
  else if (tabstyle == SPLINE) {
    fraction = (theta - tb->ang[itable]) * tb->invdelta;

    b = (theta - tb->ang[itable]) * tb->invdelta;
    a = 1.0 - b;
    u = a * tb->e[itable] + b * tb->e[itable+1] +
      ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) *
      tb->deltasq6;
    mdu_theta = a * tb->f[itable] + b * tb->f[itable+1] +
      ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
      tb->deltasq6;
    mdu_r1 = a * tb->f[itable] + b * tb->f[itable+1] +
        ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
            tb->deltasq6;
    mdu_r2 = a * tb->f[itable] + b * tb->f[itable+1] +
        ((a*a*a-a)*tb->f2[itable] + (b*b*b-b)*tb->f2[itable+1]) *
            tb->deltasq6;
  }
  */
}


/* ----------------------------------------------------------------------
   calculate potential u at angle x
------------------------------------------------------------------------- */
void AngleTable3D::u_lookup(int type, double x, double &u)
{
/*
  if (!ISFINITE(x)) {
    error->one(FLERR,"Illegal angle in angle style table");
  }

  double fraction,a,b;
  const Table *tb = &tables[tabindex[type]];
  int itable = static_cast<int> ( x * tb->invdelta);

  if (itable < 0) itable = 0;
  if (itable >= tablength) itable = tablength-1;

  if (tabstyle == LINEAR) {
    fraction = (x - tb->ang[itable]) * tb->invdelta;
    u = tb->e[itable] + fraction*tb->de[itable];
  } else if (tabstyle == SPLINE) {
    fraction = (x - tb->ang[itable]) * tb->invdelta;

    b = (x - tb->ang[itable]) * tb->invdelta;
    a = 1.0 - b;
    u = a * tb->e[itable] + b * tb->e[itable+1] +
        ((a*a*a-a)*tb->e2[itable] + (b*b*b-b)*tb->e2[itable+1]) *
            tb->deltasq6;
  }
 */
}
