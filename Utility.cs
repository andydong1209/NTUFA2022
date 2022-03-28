using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class Utils
    {
        public static void QL_Require(bool statement, string message)
        {
            if (!statement)
                throw new ApplicationException(message);
        }
    }

    public static class DStat
    {
        // The cumulative univariate normal distribution function
        public static double NormDist(double x)
        {
            // The cumulative normal distribution function            
            double z;
            if (x == 0)
                z = 0.5;
            else
            {
                double L, k;
                const double a1 = 0.31938153;
                const double a2 = -0.356563782;
                const double a3 = 1.781477937;
                const double a4 = -1.821255978;
                const double a5 = 1.330274429;

                L = Math.Abs(x);
                k = 1 / (1 + 0.2316419 * L);
                z = 1 - 1 / Math.Sqrt(2 * Math.PI) * Math.Exp(-L * L / 2)
                    * (k * (a1 + k * (a2 + k * (a3 + k * (a4 + k * a5)))));

                if (x < 0)
                    z = 1 - z;
            }
            return z;
        }

        // normal distribution function    
        public static double NormProb(double z)
        {
            return (1.0 / Math.Sqrt(2.0 * Math.PI)) * Math.Exp(-0.5 * z * z);
        }

        // The cumulative bivariate normal distribution function
        public static double NormDist(double a, double b, double rho)
        {
            double rho1, rho2, delta;
            double a1, b1, Sum;
            int i, j;
            double result = 0;

            double[] X = { 0.24840615, 0.39233107, 0.21141819, 0.03324666, 0.00082485334 };
            double[] Y = { 0.10024215, 0.48281397, 1.0609498, 1.7797294, 2.6697604 };
            a1 = a / Math.Sqrt(2 * (1 - rho * rho));
            b1 = b / Math.Sqrt(2 * (1 - rho * rho));

            if (a <= 0 && b <= 0 && rho <= 0)
            {
                Sum = 0;
                for (i = 0; i < 5; i++)
                {
                    for (j = 0; j < 5; j++)
                    {
                        Sum = Sum + X[i] * X[j] * Math.Exp(a1 * (2 * Y[i] - a1)
                            + b1 * (2 * Y[j] - b1) + 2 * rho * (Y[i] - a1) * (Y[j] - b1));
                    }
                }
                result = Math.Sqrt(1 - rho * rho) / Math.PI * Sum;
            }
            else if (a <= 0 && b >= 0 && rho >= 0)
            {
                result = NormDist(a) - NormDist(a, -b, -rho);
            }
            else if (a >= 0 && b <= 0 && rho >= 0)
            {
                result = NormDist(b) - NormDist(-a, b, -rho);
            }
            else if (a >= 0 && b >= 0 && rho <= 0)
            {
                result = NormDist(a) + NormDist(b) - 1 + NormDist(-a, -b, rho);
            }
            else if (a * b * rho > 0)
            {
                rho1 = (rho * a - b) * Math.Sign(a)
                    / Math.Sqrt(a * a - 2 * rho * a * b + b * b);
                rho2 = (rho * b - a) * Math.Sign(b)
                    / Math.Sqrt(a * a - 2 * rho * a * b + b * b);
                delta = (1 - Math.Sign(a) * Math.Sign(b)) / 4;

                result = NormDist(a, 0, rho1) + NormDist(b, 0, rho2) - delta;
            }

            return result;
        }

        public static double N_Inv(double x)
        {
            //const double SQRT_TWO_PI = 2.506628274631;
            const double e_1 = -39.6968302866538;
            const double e_2 = 220.946098424521;
            const double e_3 = -275.928510446969;
            const double e_4 = 138.357751867269;
            const double e_5 = -30.6647980661472;
            const double e_6 = 2.50662827745924;

            const double f_1 = -54.4760987982241;
            const double f_2 = 161.585836858041;
            const double f_3 = -155.698979859887;
            const double f_4 = 66.8013118877197;
            const double f_5 = -13.2806815528857;

            const double g_1 = -0.00778489400243029;
            const double g_2 = -0.322396458041136;
            const double g_3 = -2.40075827716184;
            const double g_4 = -2.54973253934373;
            const double g_5 = 4.37466414146497;
            const double g_6 = 2.93816398269878;

            const double h_1 = 0.00778469570904146;
            const double h_2 = 0.32246712907004;
            const double h_3 = 2.445134137143;
            const double h_4 = 3.75440866190742;

            const double x_l = 0.02425;
            const double x_u = 0.97575;

            double z, r;

            // Lower region: 0 < x < x_l
            if (x < x_l)
            {
                z = Math.Sqrt(-2.0 * Math.Log(x));
                z = (((((g_1 * z + g_2) * z + g_3) * z + g_4) * z + g_5) * z + g_6)
                    / ((((h_1 * z + h_2) * z + h_3) * z + h_4) * z + 1.0);
            }
            // Central region: x_l <= x <= x_u
            else if (x <= x_u)
            {
                z = x - 0.5;
                r = z * z;
                z = (((((e_1 * r + e_2) * r + e_3) * r + e_4) * r + e_5) * r + e_6) * z
                    / (((((f_1 * r + f_2) * r + f_3) * r + f_4) * r + f_5) * r + 1.0);
            }
            // Upper region. ( x_u < x < 1 )
            else
            {
                z = Math.Sqrt(-2.0 * Math.Log(1.0 - x));
                z = -(((((g_1 * z + g_2) * z + g_3) * z + g_4) * z + g_5) * z + g_6)
                    / ((((h_1 * z + h_2) * z + h_3) * z + h_4) * z + 1.0);
            }

            // Now |relative error| < 1.15e-9.  One iteration of Halley's third
            // order zero finder gives full machine precision:
            //
            //r = (N(z) - x) * SQRT_TWO_PI * exp( 0.5 * z * z );  //  f(z)/df(z)
            //z -= r/(1+0.5*z*r);

            return z;
        }

    }

    public class MersenneTwister
    {
        #region Constants -------------------------------------------------------

        private const int N = 624;
        private const int M = 397;
        private const uint MATRIX_A = 0x9908b0dfU;
        private const uint UPPER_MASK = 0x80000000U;
        private const uint LOWER_MASK = 0x7fffffffU;
        private const int MAX_RAND_INT = 0x7fffffff;

        #endregion Constants

        #region Instance Variables ----------------------------------------------

        private uint[] mag01 = { 0x0U, MATRIX_A };
        private uint[] mt = new uint[N];
        private int mti = N + 1;

        #endregion Instance Variables

        #region Constructors ----------------------------------------------------

        public MersenneTwister()
        {
            init_genrand((uint)DateTime.Now.Millisecond);
        }

        public MersenneTwister(int seed)
        {
            init_genrand((uint)seed);
        }

        public MersenneTwister(int[] init)
        {
            uint[] initArray = new uint[init.Length];
            for (int i = 0; i < init.Length; ++i)
                initArray[i] = (uint)init[i];

            init_by_array(initArray, (uint)initArray.Length);
        }

        #endregion Constructors

        #region Properties ------------------------------------------------------

        public static int MaxRandomInt
        {
            get
            {
                return 0x7fffffff;
            }
        }

        #endregion Properties

        #region Member Functions ------------------------------------------------

        public int Next()
        {
            return genrand_int31();
        }

        public int Next(int maxValue)
        {
            return Next(0, maxValue);
        }

        public int Next(int minValue, int maxValue)
        {
            if (minValue > maxValue)
            {
                int tmp = maxValue;
                maxValue = minValue;
                minValue = tmp;
            }

            return (int)(Math.Floor((maxValue - minValue + 1) * genrand_real1() + minValue));
        }

        public float NextFloat()
        {
            return (float)genrand_real2();
        }

        public float NextFloat(bool includeOne)
        {
            if (includeOne)
            {
                return (float)genrand_real1();
            }
            return (float)genrand_real2();
        }

        public float NextFloatPositive()
        {
            return (float)genrand_real3();
        }

        public double NextDouble()
        {
            return genrand_real2();
        }

        public double NextDouble(bool includeOne)
        {
            if (includeOne)
            {
                return genrand_real1();
            }
            return genrand_real2();
        }

        public double NextDoublePositive()
        {
            return genrand_real3();
        }

        public double Next53BitRes()
        {
            return genrand_res53();
        }

        public void Initialize()
        {
            init_genrand((uint)DateTime.Now.Millisecond);
        }

        public void Initialize(int seed)
        {
            init_genrand((uint)seed);
        }

        public void Initialize(int[] init)
        {
            uint[] initArray = new uint[init.Length];
            for (int i = 0; i < init.Length; ++i)
                initArray[i] = (uint)init[i];

            init_by_array(initArray, (uint)initArray.Length);
        }

        #region Methods ported from C -------------------------------------------

        private void init_genrand(uint s)
        {
            mt[0] = s & 0xffffffffU;
            for (mti = 1; mti < N; mti++)
            {
                mt[mti] =
                  (uint)(1812433253U * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
                mt[mti] &= 0xffffffffU;
            }
        }

        private void init_by_array(uint[] init_key, uint key_length)
        {
            int i, j, k;

            init_genrand(19650218U);
            i = 1; j = 0;
            k = (int)(N > key_length ? N : key_length);

            for (; k > 0; k--)
            {
                mt[i] = (uint)((uint)(mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525U)) + init_key[j] + j);
                mt[i] &= 0xffffffffU;
                i++; j++;
                if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
                if (j >= key_length) j = 0;
            }
            for (k = N - 1; k > 0; k--)
            {
                mt[i] = (uint)((uint)(mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941U)) - i);
                mt[i] &= 0xffffffffU;
                i++;
                if (i >= N) { mt[0] = mt[N - 1]; i = 1; }
            }

            mt[0] = 0x80000000U;
        }

        uint genrand_int32()
        {
            uint y;
            if (mti >= N)
            {
                int kk;

                if (mti == N + 1)
                    init_genrand(5489U);

                for (kk = 0; kk < N - M; kk++)
                {
                    y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                    mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1U];
                }
                for (; kk < N - 1; kk++)
                {
                    y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                    mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1U];
                }
                y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
                mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1U];

                mti = 0;
            }

            y = mt[mti++];
            y ^= (y >> 11);
            y ^= (y << 7) & 0x9d2c5680U;
            y ^= (y << 15) & 0xefc60000U;
            y ^= (y >> 18);

            return y;
        }

        // generates a random number on [0,0x7fffffff]-interval
        private int genrand_int31()
        {
            return (int)(genrand_int32() >> 1);
        }

        // generates a random number on [0,1]-real-interval    
        double genrand_real1()
        {
            return genrand_int32() * (1.0 / 4294967295.0);
        }

        // generates a random number on [0,1)-real-interval
        double genrand_real2()
        {
            return genrand_int32() * (1.0 / 4294967296.0);
        }

        // generates a random number on (0,1)-real-interval
        double genrand_real3()
        {
            return (((double)genrand_int32()) + 0.5) * (1.0 / 4294967296.0);
        }

        // generates a random number on [0,1) with 53-bit resolution
        double genrand_res53()
        {
            uint a = genrand_int32() >> 5, b = genrand_int32() >> 6;
            return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
        }

        #endregion Methods ported from C

        #endregion Member Functions
    }

}
