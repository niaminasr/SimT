#ifndef PAREIT_SRC_SOLVER_MATRIX
#define PAREIT_SRC_SOLVER_MATRIX

#include "ParEITConfig.hpp"
#include "geometry/Form.hpp"
#include "geometry/Grid.hpp"
#include "tools/EitBiVector.hpp"

typedef std::array<uint_t, EIT_DIM> Vertexint;
/* THIS CLASS NEEDS TO BE MODIFIED, I WILL DO SO IN THE LAST WEEK OF DECEMBER, (I have a
 * communication bug that I need to find:( )*/
class Matrix
{
   public:
    // Constructor
    Matrix( Grid *grid, Form *form, EitBiVector<real_t> *sigma ) : m_grid( grid ), m_form( form )
    {
        m_sigma = sigma;
        m_NxBeg = m_grid->getLocalBeg();
        m_NxEnd = m_grid->getLocalEnd();
    }

    // Copy Constructor
    Matrix( const Matrix &that )
        : m_grid( that.m_grid )
        , m_form( that.m_form )
        , m_sigma( that.m_sigma )
        , m_NxBeg( that.m_NxBeg )
        , m_NxEnd( that.m_NxEnd )
    {
    }

    // Destructor
    ~Matrix() {};

    Vertexint
    getIJKfromglobalidx( real_t idxglob )
    {
        Vertexint ijk;

        ijk[0] = idxglob / m_grid->getSizes()[0];
        ijk[1] = fmod( idxglob, m_grid->getSizes()[0] );

        return ijk;
    }

    /*
       The neighbors stored in class Point are determined by their global index,
       so we need to ransform them into local indices, when we work locally.
       The function Computelocalindex does so by comunicationg via a breadcast the number of nodes
       managed on each processor and then the local index is computed by substracting the effective
       of nodes on all previos processors from the global index.
    */

    real_t
    computeLocalIndex( real_t idxglob )
    {
        uint_t idxloc;
        uint_t nbPreviousPoints = 0;

        if ( EIT_PROC_RANK == EIT_MASTER_PROC )
        {
            nbPreviousPoints = 0;
        }
        else
        {
            for ( uint i = 0; i < parallel::rank; ++i )
            {
                nbPreviousPoints += m_grid->getSizeAllProc()[i];
            }
        }

        Vertexint ijk = getIJKfromglobalidx( idxglob );

        if ( ijk[0] >= m_NxBeg && ijk[0] <= m_NxEnd )
        {
            idxloc = idxglob - nbPreviousPoints;
        }
        else
        {
            idxloc = idxglob;
        }

        return idxloc;
    }

    bool
    getMyProc( real_t idx )
    {
        bool isitmine = false;

        Vertexint ijk = getIJKfromglobalidx( idx );

        if ( ( m_NxBeg <= ijk[0] ) && ( ijk[0] <= m_NxEnd ) )
        {
            isitmine = true;
        }

        return isitmine;
    }

    // Scalar product definition
    EitBiVector<real_t>
    operator*( const EitBiVector<real_t> &u )
    {
        real_t              a;
        EitBiVector<real_t> v = m_grid->createBiVector<real_t>();
        real_t              alphaj, betaj, alphai, betai, alphak, betak;
        real_t              denom;

        uint_t Nx = m_grid->getSizes()[0];
        uint_t Ny = m_grid->getSizes()[1];
        // Get the step size
        real_t hx = m_grid->getPads()[0];
        real_t hy = m_grid->getPads()[1];

        // Intialisation
        // Bradcast to all processors the number of nodes managed by each proc and stoch it into a
        // array of size ::: parallel::size. Est ce que c'est une bonne idée?

        m_grid->bcastSizes();

        MPI_Barrier( MPI_COMM_WORLD );

        const std::vector<Point *> &pts       = m_grid->getPoints();
        const std::vector<Point *> &itfs      = m_grid->getInterfaces();
        const std::vector<Point *> &haloleft  = m_grid->getHaloR();
        const std::vector<Point *> &haloright = m_grid->getHaloL();
        const uint_t               &n_pts     = pts.size();
        const uint_t               &n_itfs    = itfs.size();

        std::vector<real_t> before, after, before_sigma, after_sigma;
        after.resize( Ny );
        before.resize( Ny );
        after_sigma.resize( Ny );
        before_sigma.resize( Ny );

        // Communicate the values of "u" and "sigma" needed to write the stencil in parallel mode.
        if ( parallel::rank + 1 < parallel::size )
            MPI_Sendrecv( &u.at( size_t( n_pts - Ny ) ), Ny, EIT_MPI_REAL, parallel::rank + 1, 0,
                          &after.front(), Ny, EIT_MPI_REAL, parallel::rank + 1, 0, MPI_COMM_WORLD,
                          MPI_STATUS_IGNORE );

        if ( parallel::rank - 1 >= 0 )
            MPI_Sendrecv( &u.at( 0 ), Ny, EIT_MPI_REAL, parallel::rank - 1, 0, &before.front(), Ny,
                          EIT_MPI_REAL, parallel::rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

        if ( parallel::rank + 1 < parallel::size )
            MPI_Sendrecv( &m_sigma->at( size_t( n_pts - Ny ) ), Ny, EIT_MPI_REAL,
                          parallel::rank + 1, 0, &after_sigma.front(), Ny, EIT_MPI_REAL,
                          parallel::rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

        if ( parallel::rank - 1 >= 0 )
            MPI_Sendrecv( &m_sigma->at( 0 ), Ny, EIT_MPI_REAL, parallel::rank - 1, 0,
                          &before_sigma.front(), Ny, EIT_MPI_REAL, parallel::rank - 1, 0,
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE );

        // Initialize the stencil
        /*
        idx (or idx_itfs) is going to represent the local index of each point used in the stencil
        for the laplacien, s (or s_itfs) represents the value of the conductivity sigma at each one
        of these nodes, uij (or uij_itfs) the value of the potential "u" on ecah stencil node,
        kij(or kij_itfs)  the laplacien facteur, and finally h  the distance between the center
        stencil node uij[4] and the other ones.
        */
        uint_t  sizeStencil = 5;
        real_t  idx[sizeStencil];
        real_t  s[sizeStencil], uij[sizeStencil], kij[sizeStencil], h[sizeStencil];
        real_t *coeff = nullptr;

        real_t dx = hx * hx, dy = hy * hy;
        // i = 1 ---> Nx-2 (without the first and last column of the cartesian grid :: to set
        // boundary conditions later)
        uint_t iloc_beg = std::max( (unsigned int)1, m_NxBeg );
        uint_t iloc_end = std::min( Nx - 2, m_NxEnd );

        // The laplacien part of the matricial vector product
        for ( uint_t i = iloc_beg; i < iloc_end + 1; ++i )
        {
            for ( uint_t j = 1; j < Ny - 1; ++j )
            {
                idx[4]    = ( i - m_NxBeg ) * Nx + j;
                s[4]      = m_sigma->at( idx[4] );
                uij[4]    = u.at( idx[4] );
                real_t xi = pts[idx[4]]->getCoordinate( 0 );
                real_t yj = pts[idx[4]]->getCoordinate( 1 );

                for ( uint_t k = 0; k < sizeStencil - 1; ++k )
                {
                    idx[k] = pts[idx[4]]->getNeighs()[k];
                    if ( idx[k] > 0 )
                    {
                        /*
                         The neighbors stored in class Point are determined by their global index,
                         so we need to ransform them into local indices, when we work locally.
                         The function Computelocalindex does so by comunicationg via a breadcast the
                         number of nodes managed on each processor and then the local index is
                         computed by substracting the effective of nodes on all previos processors
                         from the global index.
                        */
                        idx[k] = computeLocalIndex( idx[k] );
                    }
                }
                // Defining Pads :: distance between the central node and the other nodes (for each
                // node). In the case of regular nodes:
                for ( uint_t k = 0; k < sizeStencil - 3; ++k )
                {
                    h[k]     = hx;
                    h[k + 2] = hy;
                }

                // Modifying pads in the case of an interface point as a neighbor
                // In the case of Irregular nodes: checked in this order West, Est, South, North.
                if ( idx[0] < 0 )
                    h[0] = xi - itfs[-idx[0] - 1]->getCoordinate( 0 );
                if ( idx[1] < 0 )
                    h[1] = itfs[-idx[1] - 1]->getCoordinate( 0 ) - xi;
                if ( idx[2] < 0 )
                    h[2] = yj - itfs[-idx[2] - 1]->getCoordinate( 1 );
                if ( idx[3] < 0 )
                    h[3] = itfs[-idx[3] - 1]->getCoordinate( 1 ) - yj;

                if ( ( i == iloc_beg ) && ( parallel::rank != EIT_MASTER_PROC ) )
                {
                    uij[0] = before.at( j );
                    s[0]   = before_sigma.at( j );

                    for ( uint_t k = 1; k < sizeStencil - 1; ++k )
                    {
                        uij[k] = u.at( idx[k] );
                        s[k]   = m_sigma->at( idx[k] );
                    }
                }
                else if ( ( i == iloc_end ) && ( parallel::rank != parallel::size - 1 ) )
                {
                    uij[0] = u.at( idx[0] );
                    s[0]   = m_sigma->at( idx[0] );
                    uij[1] = after.at( j );
                    s[1]   = after_sigma.at( j );
                    for ( uint_t k = 2; k < sizeStencil - 1; ++k )
                    {
                        uij[k] = u.at( idx[k] );
                        s[k]   = m_sigma->at( idx[k] );
                    }
                }
                else
                {
                    for ( uint_t k = 0; k < sizeStencil - 1; ++k )
                    {
                        uij[k] = u.at( idx[k] );
                        s[k]   = m_sigma->at( idx[k] );
                    }
                }

                kij[4] = 0.0;
                for ( uint_t k = 0; k < sizeStencil - 3; ++k )
                {
                    kij[k] = ( ( ( s[4] + s[k] ) / 2.0 ) / h[k] ) / ( h[0] / 2 + h[1] / 2 );
                    kij[k + 2] =
                        ( ( ( s[4] + s[k + 2] ) / 2.0 ) / h[k + 2] ) / ( h[2] / 2 + h[3] / 2 );
                    kij[4] -= ( kij[k] + kij[k + 2] );
                }

                coeff  = &v.at( idx[4] );
                *coeff = 0.;
                // for (uint_t k = 0; k < sizeStencil; ++k)
                // {
                //     *coeff -= kij[k] * uij[k];
                // }
                *coeff = 2 * u.at( idx[4] );

                // std::cout << uij[0] << "  " << uij[1] << "  " << uij[4] << "  " << uij[2] << "  "
                // << uij[3] << "  "
                //           << "i,j=" << i << "  " << j << std::endl;
                // std::cout << kij[0] << "  " << kij[1] << "  " << kij[4] << "  " << kij[2] << "  "
                // << kij[3] << std::endl;
            }
        }

        // // The Interface points lines of the matricial vector product.
        // uint_t sizeStencil_itfs = 3;
        // real_t idx_itfs[sizeStencil_itfs];
        // real_t s_itfs[sizeStencil_itfs], uij_itfs[sizeStencil_itfs], kij_itf[sizeStencil_itfs];
        // real_t x[sizeStencil_itfs], y[sizeStencil_itfs];

        // for (int_t i = 0; i < (int_t)n_itfs; ++i)
        // {

        //     idx_itfs[2] = (-(int_t)i - 1);
        //     s_itfs[2] = m_sigma->at(idx_itfs[2]);
        //     uij_itfs[2] = u.at(idx_itfs[2]);
        //     x[2] = itfs[i]->getCoordinate(0);
        //     y[2] = itfs[i]->getCoordinate(1);

        //     /* Get the stencil needed to descretizise the normal
        //     derivative on the interface pointes (the negative part of Bivector)
        //     */
        //     idx_itfs[0] = itfs[i]->getNormalStencil()[0];
        //     idx_itfs[1] = itfs[i]->getNormalStencil()[1];

        //     bool proc_p0, proc_p1;

        //     proc_p0 = getMyProc(idx_itfs[0]);
        //     proc_p1 = getMyProc(idx_itfs[1]);

        //     // If the stencil points are on the same proc
        //     if ((proc_p0 == true) && (proc_p1 == true))
        //     {

        //         idx_itfs[0] = computeLocalIndex(idx_itfs[0]);
        //         idx_itfs[1] = computeLocalIndex(idx_itfs[1]);
        //         uij_itfs[0] = u.at(idx_itfs[0]);
        //         uij_itfs[1] = u.at(idx_itfs[1]);
        //         s_itfs[0] = m_sigma->at(idx_itfs[0]);
        //         s_itfs[1] = m_sigma->at(idx_itfs[0]);

        //         // set coordinates
        //         x[0] = pts[idx_itfs[0]]->getCoordinate(0);
        //         y[0] = pts[idx_itfs[0]]->getCoordinate(1);

        //         x[1] = pts[idx_itfs[1]]->getCoordinate(0);
        //         y[1] = pts[idx_itfs[1]]->getCoordinate(1);

        //         denom = (x[2] - x[0]) * (y[2] - y[1]) - (x[2] - x[1]) * (y[2] - y[0]);
        //         alphaj = (y[0] - y[1]) / denom;
        //         betaj = (x[1] - x[0]) / denom;

        //         denom = (x[1] - x[0]) * (y[1] - y[2]) - (x[1] - x[2]) * (y[1] - y[0]);
        //         alphai = (y[0] - y[2]) / denom;
        //         betai = (x[2] - x[0]) / denom;

        //         denom = (x[0] - x[2]) * (y[0] - y[1]) - (x[0] - x[1]) * (y[0] - y[2]);
        //         alphak = (y[2] - y[1]) / denom;
        //         betak = (x[1] - x[2]) / denom;

        //         kij_itf[0] = s_itfs[0] * (alphai * itfs[i]->getNormal()[0] + betai *
        //         itfs[i]->getNormal()[1]); kij_itf[1] = s_itfs[1] * (alphak *
        //         itfs[i]->getNormal()[0] + betak * itfs[i]->getNormal()[1]); kij_itf[2] =
        //         s_itfs[2] * (alphaj * itfs[i]->getNormal()[0] + betaj * itfs[i]->getNormal()[1]);

        //         coeff = &v.at(idx_itfs[2]);
        //         *coeff = 1 * u.at(idx_itfs[2]);
        //         // if (itfs[-idx_itfs[2] - 1]->getisOnElectrode() == -1)
        //         // {
        //         //     *coeff = 0.;
        //         //     for (uint_t k = 0; k < sizeStencil_itfs; ++k)
        //         //         *coeff += kij_itf[k] * uij_itfs[k];
        //         // }
        //         // else
        //         // {
        //         //     kij_itf[2] += 1;

        //         //     *coeff = 0.;
        //         //     for (uint_t k = 0; k < sizeStencil_itfs; ++k)
        //         //         *coeff += kij_itf[k] * uij_itfs[k];

        //         //     *coeff += -1 * u.at(-(n_itfs + itfs[-idx_itfs[2] - 1]->getisOnElectrode()
        //         + 1));
        //         // }
        //     }
        //     // I am on Proc p0
        //     else if ((proc_p0 == true) && (proc_p1 == false))
        //     {
        //         // First stencil Node :: we already have it
        //         idx_itfs[0] = computeLocalIndex(idx_itfs[0]);
        //         uij_itfs[0] = u.at(idx_itfs[0]);
        //         x[0] = pts[idx_itfs[0]]->getCoordinate(0);
        //         y[0] = pts[idx_itfs[0]]->getCoordinate(1);
        //         s_itfs[0] = m_sigma->at(idx_itfs[0]);

        //         // Second stencil Node :: we need it
        //         // Here i need the index we need in the communicated vector
        //         idx_itfs[1] = getIJKfromglobalidx(idx_itfs[1])[1];

        //         // In this case I have to receive "u" and "sigma" from the next Processor
        //         if (itfs[i]->getNormalStencil()[0] < itfs[i]->getNormalStencil()[1])
        //         {
        //             uij_itfs[1] = after.at(idx_itfs[1]);
        //             s_itfs[1] = after_sigma.at(idx_itfs[1]);

        //             x[1] = haloleft[idx_itfs[1]]->getCoordinate(0);
        //             y[1] = haloleft[idx_itfs[1]]->getCoordinate(1);
        //         }
        //         // In this case I have to receive "u" and "sigma" from the previous Processor
        //         else
        //         {
        //             uij_itfs[1] = before.at(idx_itfs[1]);
        //             s_itfs[1] = before_sigma.at(idx_itfs[1]);

        //             x[1] = haloright[idx_itfs[1]]->getCoordinate(0);
        //             y[1] = haloright[idx_itfs[1]]->getCoordinate(1);
        //         }

        //         denom = (x[2] - x[0]) * (y[2] - y[1]) - (x[2] - x[1]) * (y[2] - y[0]);
        //         alphaj = (y[0] - y[1]) / denom;
        //         betaj = (x[1] - x[0]) / denom;

        //         denom = (x[1] - x[0]) * (y[1] - y[2]) - (x[1] - x[2]) * (y[1] - y[0]);
        //         alphai = (y[0] - y[2]) / denom;
        //         betai = (x[2] - x[0]) / denom;

        //         denom = (x[0] - x[2]) * (y[0] - y[1]) - (x[0] - x[1]) * (y[0] - y[2]);
        //         alphak = (y[2] - y[1]) / denom;
        //         betak = (x[1] - x[2]) / denom;

        //         kij_itf[0] = s_itfs[0] * (alphai * itfs[i]->getNormal()[0] + betai *
        //         itfs[i]->getNormal()[1]); kij_itf[1] = s_itfs[1] * (alphak *
        //         itfs[i]->getNormal()[0] + betak * itfs[i]->getNormal()[1]); kij_itf[2] =
        //         s_itfs[2] * (alphaj * itfs[i]->getNormal()[0] + betaj * itfs[i]->getNormal()[1]);

        //         coeff = &v.at(idx_itfs[2]);
        //         // if (itfs[i]->getisOnElectrode() == -1)
        //         // {
        //         //     *coeff = 0.;
        //         //     for (uint_t k = 0; k < sizeStencil_itfs; ++k)
        //         //         *coeff += kij_itf[k] * uij_itfs[k];
        //         // }
        //         // else
        //         // {
        //         //     kij_itf[2] += 1;

        //         //     *coeff = 0.;
        //         //     for (uint_t k = 0; k < sizeStencil; ++k)
        //         //         *coeff += kij_itf[k] * uij_itfs[k];

        //         //     *coeff += -1 * u.at(-(n_itfs + itfs[-idx_itfs[2] - 1]->getisOnElectrode()
        //         + 1));
        //         // }
        //         *coeff = 2 * u.at(idx_itfs[2]);
        //     }
        //     // With the help of MPI_Allreduce everyone is going to have the values of "v"
        //     calculated on each processor. MPI_Allreduce(&v.at(idx_itfs[2]), &v.at(idx_itfs[2]),
        //     1, EIT_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        // }

        // // The Electrode lines of the matricial vector product
        // uint_t nb_elec = EIT_ElEC;
        // real_t idx_elec_0, uij_elec_0, kij_elec_0;
        // real_t uij_elec[nb_elec][sizeStencil - 1], kij_elec[nb_elec][sizeStencil - 1];
        // real_t Um[nb_elec] = {0, 0, 0, 0};

        // for (uint_t k = 0; k < nb_elec; ++k)
        // {
        //     for (uint_t l = 0; l < 4; ++l)
        //     {
        //         uij_elec[k][l] = 0;
        //         kij_elec[k][l] = 0;
        //     }
        // }

        // for (uint_t i = iloc_beg; i < iloc_end + 1; ++i)
        // {

        //     for (uint_t j = 1; j < Ny - 1; ++j)
        //     {

        //         idx_elec_0 = ((i)-m_NxBeg) * Nx + j; // local index

        //         if (pts[idx_elec_0]->isIrregular() == true) // If I am an irregular point.
        //         {

        //             for (uint_t k = 0; k < sizeStencil - 1; ++k) // check the four neighbors of
        //             the point object for interface points
        //             {
        //                 if (pts[idx_elec_0]->getNeighs()[k] < 0) //  If I have interface points
        //                 as a neighbor.
        //                 {
        //                     if (itfs[-(pts[idx_elec_0]->getNeighs()[k]) - 1]->getisOnElectrode()
        //                     != -1) // If these interface points are on an electrode.
        //                     {
        //                         uint_t electrode_number = itfs[-pts[idx_elec_0]->getNeighs()[k] -
        //                         1]->getisOnElectrode(); uij_elec[electrode_number][k] =
        //                         u.at(idx_elec_0); kij_elec[electrode_number][k] = -hx * hy *
        //                         m_form->approximateIntegral(*pts[idx_elec_0], hy);
        //                         Um[electrode_number] += hx * hy *
        //                         m_form->approximateIntegral(*pts[idx_elec_0], hy);
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }

        // for (uint_t i = 0; i < nb_elec; i++)
        // {
        //     coeff = &v.at(-(n_itfs + i + 1));
        //     *coeff = 1 * u.at(-(n_itfs + i + 1));

        //     // uij_elec_0 = u.at(-(n_itfs + i + 1));
        //     // kij_elec_0 = Um[i];

        //     // *coeff = uij_elec_0 * kij_elec_0;

        //     // for (uint_t k = 0; k < sizeStencil - 1; ++k)
        //     //     *coeff += kij_elec[i][k] * uij_elec[i][k];
        //     // if (i == 0)
        //     // {
        //     //     *coeff += 0.0000000001 * u.at(-(n_itfs + 0 + 1));
        //     // }
        // }

        // Boundary Points
        // Right boundary of the cube.
        uint_t sizeStencil_bound = 4;
        real_t idx_bound[sizeStencil_bound];
        real_t s_bound[sizeStencil_bound], uij_bound[sizeStencil_bound],
            kij_bound[sizeStencil_bound], h_bound[sizeStencil_bound];

        if ( parallel::rank == EIT_MASTER_PROC )
        {
            for ( uint j = 1; j < Ny - 1; ++j )
            {
                idx_bound[3] = ( 0 - m_NxBeg ) * Nx + j;
                s_bound[3]   = m_sigma->at( idx_bound[3] );
                uij_bound[3] = u.at( idx_bound[3] );

                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    idx_bound[k] = pts[idx_bound[3]]->getNeighs()[k + 1];
                    if ( idx_bound[k] > 0 )
                    {
                        /*
                         The neighbors stored in class Point are determined by their global index,
                         so we need to ransform them into local indices, when we work locally.
                         The function Computelocalindex does so by comunicationg via a breadcast the
                         number of nodes managed on each processor and then the local index is
                         computed by substracting the effective of nodes on all previos processors
                         from the global index.
                        */
                        idx_bound[k] = computeLocalIndex( idx_bound[k] );
                    }
                }
                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    uij_bound[k] = u.at( idx_bound[k] );
                    s_bound[k]   = m_sigma->at( idx_bound[k] );
                }
                h_bound[0] = hx;
                for ( uint_t k = 1; k < sizeStencil - 1; ++k )
                {
                    h_bound[k] = hy;
                }

                kij_bound[3] = 0.0;
                kij_bound[0] =
                    ( ( ( s_bound[3] + s_bound[0] ) / 2.0 ) / h_bound[0] ) / ( h_bound[0] );
                kij_bound[1] = ( ( ( s_bound[3] + s_bound[1] ) / 2.0 ) / h_bound[1] ) /
                               ( h_bound[1] / 2 + h_bound[2] / 2 );
                kij_bound[2] = ( ( ( s_bound[3] + s_bound[2] ) / 2.0 ) / h_bound[2] ) /
                               ( h_bound[1] / 2 + h_bound[2] / 2 );

                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    kij_bound[3] -= ( kij_bound[k] );
                }
                kij_bound[3] += ( -1. / ( hy * ( 1. ) ) ) / hx;

                coeff = &v.at( idx_bound[3] );
                // *coeff = 0.;
                // for (uint_t k = 0; k < sizeStencil_bound; ++k)
                // {
                //     *coeff -= kij_bound[k] * uij_bound[k];
                // }
                *coeff = u.at( idx_bound[3] );
            }
        }

        // Left boundary of the cube
        if ( parallel::rank == parallel::size - 1 )
        {
            for ( uint j = 1; j < Ny - 1; ++j )
            {
                idx_bound[3] = ( m_NxEnd - m_NxBeg ) * Nx + j;
                s_bound[3]   = m_sigma->at( idx_bound[3] );
                uij_bound[3] = u.at( idx_bound[3] );

                idx_bound[0] = pts[idx_bound[3]]->getNeighs()[0];
                for ( uint_t k = 1; k < sizeStencil_bound - 1; ++k )
                {
                    idx_bound[k] = pts[idx_bound[3]]->getNeighs()[k + 1];
                }
                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    if ( idx_bound[k] > 0 )
                    {
                        /*
                         The neighbors stored in class Point are determined by their global index,
                         so we need to ransform them into local indices, when we work locally.
                         The function Computelocalindex does so by comunicationg via a breadcast the
                         number of nodes managed on each processor and then the local index is
                         computed by substracting the effective of nodes on all previos processors
                         from the global index.
                        */
                        idx_bound[k] = computeLocalIndex( idx_bound[k] );
                    }
                }
                // std::cout << idx_bound[0] << "  " << idx_bound[1] << "  " << idx_bound[2] << "  "
                // << idx_bound[3] << "  " << std::endl;

                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    uij_bound[k] = u.at( idx_bound[k] );
                    s_bound[k]   = m_sigma->at( idx_bound[k] );
                }
                h_bound[0] = hx;
                for ( uint_t k = 1; k < sizeStencil - 1; ++k )
                {
                    h_bound[k] = hy;
                }

                kij_bound[3] = 0.0;
                kij_bound[0] =
                    ( ( ( s_bound[3] + s_bound[0] ) / 2.0 ) / h_bound[0] ) / ( h_bound[0] );
                kij_bound[1] = ( ( ( s_bound[3] + s_bound[1] ) / 2.0 ) / h_bound[1] ) /
                               ( h_bound[1] / 2 + h_bound[2] / 2 );
                kij_bound[2] = ( ( ( s_bound[3] + s_bound[2] ) / 2.0 ) / h_bound[2] ) /
                               ( h_bound[1] / 2 + h_bound[2] / 2 );

                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    kij_bound[3] -= ( kij_bound[k] );
                }
                kij_bound[3] += ( -1. / ( hy * ( 1. ) ) ) / hx;

                coeff  = &v.at( idx_bound[3] );
                *coeff = 0.;
                // for (uint_t k = 0; k < sizeStencil_bound; ++k)
                // {
                //     *coeff -= kij_bound[k] * uij_bound[k];
                // }
                // std::cout << kij_bound[0] << "  " << kij_bound[1] << "  " << kij_bound[2] << "  "
                // << kij_bound[3] << "  " << std::endl; std::cout << uij_bound[0] << "  " <<
                // uij_bound[1] << "  " << uij_bound[2] << "  " << uij_bound[3] << "  " <<
                // std::endl;
                *coeff = u.at( idx_bound[3] );
            }
        }

        // Lower boundary of the cube
        for ( uint i = iloc_beg; i < iloc_end + 1; ++i )
        {
            idx_bound[3] = ( i - m_NxBeg ) * Nx + 0;
            s_bound[3]   = m_sigma->at( idx_bound[3] );
            uij_bound[3] = u.at( idx_bound[3] );

            idx_bound[2] = pts[idx_bound[3]]->getNeighs()[3];
            for ( uint_t k = 0; k < sizeStencil_bound - 2; ++k )
            {
                idx_bound[k] = pts[idx_bound[3]]->getNeighs()[k];
            }
            for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
            {
                if ( idx_bound[k] > 0 )
                {
                    /*
                     The neighbors stored in class Point are determined by their global index,
                     so we need to ransform them into local indices, when we work locally.
                     The function Computelocalindex does so by comunicationg via a breadcast the
                     number of nodes managed on each processor and then the local index is computed
                     by substracting the effective of nodes on all previos processors from the
                     global index.
                    */
                    idx_bound[k] = computeLocalIndex( idx_bound[k] );
                }
            }
            // std::cout << idx_bound[0] << "  " << idx_bound[1] << "  " << idx_bound[2] << "  " <<
            // idx_bound[3] << "  " << parallel::rank << std::endl;
            if ( ( i == iloc_beg ) && ( parallel::rank != EIT_MASTER_PROC ) )
            {
                uij_bound[0] = before.at( 0 );
                s_bound[0]   = before_sigma.at( 0 );

                for ( uint_t k = 1; k < sizeStencil_bound - 1; ++k )
                {
                    uij_bound[k] = u.at( idx_bound[k] );
                    s_bound[k]   = m_sigma->at( idx_bound[k] );
                }
            }
            else if ( ( i == iloc_end ) && ( parallel::rank != parallel::size - 1 ) )
            {
                uij_bound[0] = u.at( idx_bound[0] );
                s_bound[0]   = m_sigma->at( idx_bound[0] );
                uij_bound[1] = after.at( 0 );
                s_bound[1]   = after_sigma.at( 0 );
                uij_bound[2] = u.at( idx_bound[2] );
                s_bound[2]   = m_sigma->at( idx_bound[2] );
            }
            else
            {
                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    uij_bound[k] = u.at( idx_bound[k] );
                    s_bound[k]   = m_sigma->at( idx_bound[k] );
                }
            }

            h_bound[0] = hy;
            for ( uint_t k = 0; k < sizeStencil - 2; ++k )
            {
                h_bound[k] = hx;
            }

            kij_bound[3] = 0.0;
            kij_bound[0] = ( ( ( s_bound[3] + s_bound[0] ) / 2.0 ) / h_bound[0] ) /
                           ( h_bound[0] / 2 + h_bound[1] / 2 );
            kij_bound[1] = ( ( ( s_bound[3] + s_bound[1] ) / 2.0 ) / h_bound[1] ) /
                           ( h_bound[0] / 2 + h_bound[1] / 2 );
            kij_bound[2] = ( ( ( s_bound[3] + s_bound[2] ) / 2.0 ) / h_bound[2] ) / ( h_bound[2] );

            for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
            {
                kij_bound[3] -= ( kij_bound[k] );
            }
            kij_bound[3] += -1. / ( hx * hx );

            coeff  = &v.at( idx_bound[3] );
            *coeff = 0.;
            // for (uint_t k = 0; k < sizeStencil_bound; ++k)
            // {
            //     *coeff -= kij_bound[k] * uij_bound[k];
            // }
            // std::cout << uij_bound[0] << "  " << uij_bound[1] << "  " << uij_bound[2] << "  " <<
            // uij_bound[3] << "  " << parallel::rank << std::endl; std::cout << kij_bound[0] << "
            // " << kij_bound[1] << "  " << kij_bound[2] << "  " << kij_bound[3] << "  " <<
            // parallel::rank << std::endl;
            *coeff = u.at( idx_bound[3] );
        }

        // upper boundary of the cube
        for ( uint i = iloc_beg; i < iloc_end + 1; ++i )
        {
            idx_bound[3] = ( i - m_NxBeg ) * Nx + ( Ny - 1 );
            s_bound[3]   = m_sigma->at( idx_bound[3] );
            uij_bound[3] = u.at( idx_bound[3] );

            idx_bound[2] = pts[idx_bound[3]]->getNeighs()[2];
            for ( uint_t k = 0; k < sizeStencil_bound - 2; ++k )
            {
                idx_bound[k] = pts[idx_bound[3]]->getNeighs()[k];
            }
            for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
            {
                if ( idx_bound[k] > 0 )
                {
                    /*
                     The neighbors stored in class Point are determined by their global index,
                     so we need to ransform them into local indices, when we work locally.
                     The function Computelocalindex does so by comunicationg via a breadcast the
                     number of nodes managed on each processor and then the local index is computed
                     by substracting the effective of nodes on all previos processors from the
                     global index.
                    */
                    idx_bound[k] = computeLocalIndex( idx_bound[k] );
                }
            }
            // std::cout << idx_bound[0] << "  " << idx_bound[1] << "  " << idx_bound[2] << "  " <<
            // idx_bound[3] << "  " << parallel::rank<< std::endl;

            if ( ( i == iloc_beg ) && ( parallel::rank != EIT_MASTER_PROC ) )
            {
                uij_bound[0] = before.at( ( Ny - 1 ) );
                s_bound[0]   = before_sigma.at( ( Ny - 1 ) );

                for ( uint_t k = 1; k < sizeStencil_bound - 1; ++k )
                {
                    uij_bound[k] = u.at( idx_bound[k] );
                    s_bound[k]   = m_sigma->at( idx_bound[k] );
                }
            }
            else if ( ( i == iloc_end ) && ( parallel::rank != parallel::size - 1 ) )
            {
                uij_bound[0] = u.at( idx_bound[0] );
                s_bound[0]   = m_sigma->at( idx_bound[0] );
                uij_bound[1] = after.at( ( Ny - 1 ) );
                s_bound[1]   = after_sigma.at( ( Ny - 1 ) );
                uij_bound[2] = u.at( idx_bound[2] );
                s_bound[2]   = m_sigma->at( idx_bound[2] );
            }
            else
            {
                for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
                {
                    uij_bound[k] = u.at( idx_bound[k] );
                    s_bound[k]   = m_sigma->at( idx_bound[k] );
                }
            }
            h_bound[0] = hy;
            for ( uint_t k = 0; k < sizeStencil - 2; ++k )
            {
                h_bound[k] = hx;
            }

            kij_bound[3] = 0.0;
            kij_bound[0] = ( ( ( s_bound[3] + s_bound[0] ) / 2.0 ) / h_bound[0] ) /
                           ( h_bound[0] / 2 + h_bound[1] / 2 );
            kij_bound[1] = ( ( ( s_bound[3] + s_bound[1] ) / 2.0 ) / h_bound[1] ) /
                           ( h_bound[0] / 2 + h_bound[1] / 2 );
            kij_bound[2] = ( ( ( s_bound[3] + s_bound[2] ) / 2.0 ) / h_bound[2] ) / ( h_bound[2] );

            for ( uint_t k = 0; k < sizeStencil_bound - 1; ++k )
            {
                kij_bound[3] -= ( kij_bound[k] );
            }
            kij_bound[3] += -1. / ( hy * hx );

            coeff  = &v.at( idx_bound[3] );
            *coeff = 0.;
            // for (uint_t k = 0; k < sizeStencil_bound; ++k)
            // {
            //     *coeff -= kij_bound[k] * uij_bound[k];
            // }
            // std::cout << uij_bound[0] << "  " << uij_bound[1] << "  " << uij_bound[2] << "  " <<
            // uij_bound[3] << "  " << parallel::rank << std::endl; std::cout << kij_bound[0] << "
            // " << kij_bound[1] << "  " << kij_bound[2] << "  " << kij_bound[3] << "  " <<
            // parallel::rank << std::endl;
            *coeff = u.at( idx_bound[3] );
        }
        // setting solution for Point :: i=0,j=0 to 1
        coeff = &v.at( 0 );
        //*coeff = (1. / (hy * hx)) * u.at(1) + (1. / (hy * hx)) * u.at(Nx) - (4 * (1. / (hy * hx)))
        //* u.at(0);
        *coeff = 1 * u.at( 0 );
        // Point :: i=0,j=Ny
        coeff = &v.at( Ny - 1 );
        //*coeff = (1. / (hy * hx)) * u.at(Ny-2) + (1. / (hy * hx)) * u.at(Ny-1+Nx) - (4 * (1. / (hy
        //* hx))) * u.at(Ny-1);
        *coeff = u.at( Ny - 1 );
        // Point :: i=Nx,j=0
        coeff = &v.at( ( Nx - 1 ) * Nx );
        //*coeff = (1. / (hy * hx)) * u.at((Nx-1)*Nx-Nx) + (1. / (hy * hx)) * u.at((Nx-1)*Nx + 1) -
        //(4 * (1. / (hy * hx))) * u.at((Nx-1)*Nx);
        *coeff = u.at( ( Nx - 1 ) * Nx );
        // Point ::i=Nx,j=Ny
        coeff = &v.at( ( Nx - 1 ) * Nx + Ny - 1 );
        //*coeff = (1. / (hy * hx)) * u.at((Nx-1)*Nx +Ny-1 - Nx) + (1. / (hy * hx)) * u.at((Nx-1)*Nx
        //+Ny-1 - 1) - (4 * (1. / (hy * hx))) * u.at((Nx-1)*Nx +Ny-1);
        *coeff = u.at( ( Nx - 1 ) * Nx + Ny - 1 );

        return v;
    }

   private:
    Grid                *m_grid;
    Form                *m_form;
    EitBiVector<real_t> *m_sigma;

    uint_t m_NxBeg, m_NxEnd;
};

// EitBiVector<real_t> operator* (const Matrix &A, const EitBiVector<real_t> &u);

#endif /* PAREIT_SRC_SOLVER_MATRIX */
