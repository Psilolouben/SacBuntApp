using System;
using System.Collections.Generic;

namespace KlarityRisk.ExecuteJob.EvaluationMethods
{
    internal static class MathHelpers
    {

        public static double LinearInterpolation(List<double> knownSamples, double k)
        {
            int k_hat = (int)Math.Floor(k);
            double a = k - k_hat;
            //double f_k = a * knownSamples[k_hat + 1] + (1 - a) * knownSamples[k_hat];
            double f_k = a * knownSamples[k_hat] + (1 - a) * knownSamples[k_hat - 1];
            return f_k;
        }

        public static Tuple<int, double> LinearInterpolation(List<Tuple<int, double>> knownSamples, double k)
        {
            int k_hat = (int)Math.Floor(k);
            double a = k - k_hat;
            //double f_k = a * knownSamples[k_hat + 1] + (1 - a) * knownSamples[k_hat];
            Tuple<int, double> f_k = new Tuple<int, double>(knownSamples[k_hat].Item1, a * knownSamples[k_hat].Item2 + (1 - a) * knownSamples[k_hat - 1].Item2);
            return f_k;
        }

        internal static void GetHistogramValues(double[] R, int N, double k, ref double[] Hist, ref double[] Intervals, ref int k_star)
        {
            int b = (int)Math.Min(Math.Round((decimal)N / 15), 50);

            Intervals = new double[b + 1];

            for (int l = 0; l <= b; l++)
            {
                Intervals[l] = R[0] + (R[N - 1] - R[0]) * ((double)l / b);

            }

            Intervals[b] = R[R.Length - 1];

            Hist = new double[b];
            int i = 0; int j = 0;
            bool Allocated = false;
            int k_round = (int)Math.Round(k);

            while (i < N)
            {
                Allocated = false;

                while (Allocated == false)
                {
                    if (R[i] <= Intervals[j + 1])
                    {
                        Hist[j] = Hist[j] + 1;
                        Allocated = true;
                        if (i == k_round)
                        {
                            k_star = j;
                        }
                    }
                    else
                    {
                        j++;
                    }
                }
                i++;
            }
            k_star++;
        }

        public class SimpleRNG
        {
            private static uint m_w;
            private static uint m_z;

            static SimpleRNG()
            {
                // These values are not magical, just the default values Marsaglia used.
                // Any pair of unsigned integers should be fine.
                m_w = 521288629;
                m_z = 362436069;
            }

            // The random generator seed can be set three ways:
            // 1) specifying two non-zero unsigned integers
            // 2) specifying one non-zero unsigned integer and taking a default value for the second
            // 3) setting the seed from the system time

            public static void SetSeed(uint u, uint v)
            {
                if (u != 0) m_w = u;
                if (v != 0) m_z = v;
            }

            public static void SetSeed(uint u)
            {
                m_w = u;
            }

            public static void SetSeedFromSystemTime()
            {
                System.DateTime dt = System.DateTime.Now;
                long x = dt.ToFileTime();
                SetSeed((uint)(x >> 16), (uint)(x % 4294967296));
            }

            // Produce a uniform random sample from the open interval (0, 1).
            // The method will not return either end point.
            public static double GetUniform()
            {
                // 0 <= u < 2^32
                uint u = GetUint();
                // The magic number below is 1/(2^32 + 2).
                // The result is strictly between 0 and 1.
                return (u + 1.0) * 2.328306435454494e-10;
            }

            // This is the heart of the generator.
            // It uses George Marsaglia's MWC algorithm to produce an unsigned integer.
            // See http://www.bobwheeler.com/statistics/Password/MarsagliaPost.txt
            private static uint GetUint()
            {
                m_z = 36969 * (m_z & 65535) + (m_z >> 16);
                m_w = 18000 * (m_w & 65535) + (m_w >> 16);
                return (m_z << 16) + m_w;
            }

            // Get normal (Gaussian) random sample with mean 0 and standard deviation 1
            public static double GetNormal()
            {
                // Use Box-Muller algorithm
                double u1 = GetUniform();
                double u2 = GetUniform();
                double r = Math.Sqrt(-2.0 * Math.Log(u1));
                double theta = 2.0 * Math.PI * u2;
                return r * Math.Sin(theta);
            }

            // Get normal (Gaussian) random sample with specified mean and standard deviation
            public static double GetNormal(double mean, double standardDeviation)
            {
                if (standardDeviation <= 0.0)
                {
                    string msg = string.Format("Shape must be positive. Received {0}.", standardDeviation);
                    throw new ArgumentOutOfRangeException(msg);
                }
                return mean + standardDeviation * GetNormal();
            }

            // Get exponential random sample with mean 1
            public static double GetExponential()
            {
                return -Math.Log(GetUniform());
            }

            // Get exponential random sample with specified mean
            public static double GetExponential(double mean)
            {
                if (mean <= 0.0)
                {
                    string msg = string.Format("Mean must be positive. Received {0}.", mean);
                    throw new ArgumentOutOfRangeException(msg);
                }
                return mean * GetExponential();
            }

            public static double GetGamma(double shape, double scale)
            {
                // Implementation based on "A Simple Method for Generating Gamma Variables"
                // by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
                // Vol 26, No 3, September 2000, pages 363-372.

                double d, c, x, xsquared, v, u;

                if (shape >= 1.0)
                {
                    d = shape - 1.0 / 3.0;
                    c = 1.0 / Math.Sqrt(9.0 * d);
                    for (;;)
                    {
                        do
                        {
                            x = GetNormal();
                            v = 1.0 + c * x;
                        }
                        while (v <= 0.0);
                        v = v * v * v;
                        u = GetUniform();
                        xsquared = x * x;
                        if (u < 1.0 - .0331 * xsquared * xsquared || Math.Log(u) < 0.5 * xsquared + d * (1.0 - v + Math.Log(v)))
                            return scale * d * v;
                    }
                }
                else if (shape <= 0.0)
                {
                    string msg = string.Format("Shape must be positive. Received {0}.", shape);
                    throw new ArgumentOutOfRangeException(msg);
                }
                else
                {
                    double g = GetGamma(shape + 1.0, 1.0);
                    double w = GetUniform();
                    return scale * g * Math.Pow(w, 1.0 / shape);
                }
            }


        }

        internal static double COV(int x, int y, double[,] Returns, int t, bool useExponentialVolatility, double lamda)
        {
            Double Sum = 0.0;
            Double CurrentCOV = 0.0;
            var nonInfiniteValues = default(double);


            if (!useExponentialVolatility) //use equation 13
            {
                for (int i = 0; i < t; i++)
                {
                    if (!double.IsNaN(Returns[i, x]) && !double.IsNaN(Returns[i, y]) && !double.IsInfinity(Returns[i, x]) && !double.IsInfinity(Returns[i, y]))
                    {
                        Sum += Returns[i, x] * Returns[i, y];
                        nonInfiniteValues++;

                    }

                }


                CurrentCOV = (1.0 / (nonInfiniteValues)) * Sum;
            }
            else
            {
                double inner_sum = 0;
                for (int i = 0; i < t; i++)
                {
                    if (!double.IsNaN(Returns[i, x]) && !double.IsNaN(Returns[i, y]) && !double.IsInfinity(Returns[i, x]) && !double.IsInfinity(Returns[i, y]))
                    {
                        inner_sum += Math.Pow(lamda, i) * Returns[i, x] * Returns[i, y];
                        nonInfiniteValues++;
                    }
                }
                CurrentCOV = (1 - lamda) * inner_sum;
            }

            return CurrentCOV;
        }

        internal static double COV_Stress(int x, int y, double[,] Returns, int t, bool useExponentialVolatility, double lamda)
        {
            Double Sum = 0.0;

            Double CurrentCOV = 0.0;
            var nonInfiniteValues = default(double);

            if (!useExponentialVolatility) //use equation 13
            {
                for (int i = 0; i < t; i++)
                {
                    if (!double.IsNaN(Returns[i, x]) && !double.IsNaN(Returns[i, y]) && !double.IsInfinity(Returns[i, x]) && !double.IsInfinity(Returns[i, y]))
                    {
                        Sum += Returns[i, x] * Returns[i, y];
                        nonInfiniteValues++;
                    }
                }

                CurrentCOV = (double)Sum / t;
            }
            else
            {
                double inner_sum = 0;
                for (int i = 0; i < t; i++)
                {
                    if (!double.IsNaN(Returns[i, x]) && !double.IsNaN(Returns[i, y]) && !double.IsInfinity(Returns[i, x]) && !double.IsInfinity(Returns[i, y]))
                    {
                        inner_sum += Math.Pow(lamda, i) * Returns[i, x] * Returns[i, y];
                        nonInfiniteValues++;
                    }
                }
                CurrentCOV = (1 - lamda) * inner_sum;
            }

            return CurrentCOV;
        }



        internal static double CND(double X)
        {
            double L = 0.0;
            double K = 0.0;
            double dCND = 0.0;
            const double a1 = 0.31938153;
            const double a2 = -0.356563782;
            const double a3 = 1.781477937;
            const double a4 = -1.821255978;
            const double a5 = 1.330274429;
            L = Math.Abs(X);
            K = 1.0 / (1.0 + 0.2316419 * L);
            dCND = 1.0 - 1.0 / Math.Sqrt(2 * Convert.ToDouble(Math.PI.ToString())) *
                Math.Exp(-L * L / 2.0) * (a1 * K + a2 * K * K + a3 * Math.Pow(K, 3.0) +
                a4 * Math.Pow(K, 4.0) + a5 * Math.Pow(K, 5.0));

            if (X < 0)
            {
                return 1.0 - dCND;
            }
            else
            {
                return dCND;
            }
        }

        internal static double inverf(double e)
        {
            double result = 0;

            result = invnormaldistribution(0.5 * (e + 1)) / Math.Sqrt(2);
            return result;
        }

        internal static double invnormaldistribution(double y0)
        {
            double result = 0;
            double expm2 = 0;
            double s2pi = 0;
            double x = 0;
            double y = 0;
            double z = 0;
            double y2 = 0;
            double x0 = 0;
            double x1 = 0;
            int code = 0;
            double p0 = 0;
            double q0 = 0;
            double p1 = 0;
            double q1 = 0;
            double p2 = 0;
            double q2 = 0;

            expm2 = 0.13533528323661269189;
            s2pi = 2.50662827463100050242;
            if ((double)(y0) <= (double)(0))
            {
                result = double.NegativeInfinity;
                return result;
            }
            if ((double)(y0) >= (double)(1))
            {
                result = double.PositiveInfinity;
                return result;
            }
            code = 1;
            y = y0;
            if ((double)(y) > (double)(1.0 - expm2))
            {
                y = 1.0 - y;
                code = 0;
            }
            if ((double)(y) > (double)(expm2))
            {
                y = y - 0.5;
                y2 = y * y;
                p0 = -59.9633501014107895267;
                p0 = 98.0010754185999661536 + y2 * p0;
                p0 = -56.6762857469070293439 + y2 * p0;
                p0 = 13.9312609387279679503 + y2 * p0;
                p0 = -1.23916583867381258016 + y2 * p0;
                q0 = 1;
                q0 = 1.95448858338141759834 + y2 * q0;
                q0 = 4.67627912898881538453 + y2 * q0;
                q0 = 86.3602421390890590575 + y2 * q0;
                q0 = -225.462687854119370527 + y2 * q0;
                q0 = 200.260212380060660359 + y2 * q0;
                q0 = -82.0372256168333339912 + y2 * q0;
                q0 = 15.9056225126211695515 + y2 * q0;
                q0 = -1.18331621121330003142 + y2 * q0;
                x = y + y * y2 * p0 / q0;
                x = x * s2pi;
                result = x;
                return result;
            }
            x = Math.Sqrt(-(2.0 * Math.Log(y)));
            x0 = x - Math.Log(x) / x;
            z = 1.0 / x;
            if ((double)(x) < (double)(8.0))
            {
                p1 = 4.05544892305962419923;
                p1 = 31.5251094599893866154 + z * p1;
                p1 = 57.1628192246421288162 + z * p1;
                p1 = 44.0805073893200834700 + z * p1;
                p1 = 14.6849561928858024014 + z * p1;
                p1 = 2.18663306850790267539 + z * p1;
                p1 = -(1.40256079171354495875 * 0.1) + z * p1;
                p1 = -(3.50424626827848203418 * 0.01) + z * p1;
                p1 = -(8.57456785154685413611 * 0.0001) + z * p1;
                q1 = 1;
                q1 = 15.7799883256466749731 + z * q1;
                q1 = 45.3907635128879210584 + z * q1;
                q1 = 41.3172038254672030440 + z * q1;
                q1 = 15.0425385692907503408 + z * q1;
                q1 = 2.50464946208309415979 + z * q1;
                q1 = -(1.42182922854787788574 * 0.1) + z * q1;
                q1 = -(3.80806407691578277194 * 0.01) + z * q1;
                q1 = -(9.33259480895457427372 * 0.0001) + z * q1;
                x1 = z * p1 / q1;
            }
            else
            {
                p2 = 3.23774891776946035970;
                p2 = 6.91522889068984211695 + z * p2;
                p2 = 3.93881025292474443415 + z * p2;
                p2 = 1.33303460815807542389 + z * p2;
                p2 = 2.01485389549179081538 * 0.1 + z * p2;
                p2 = 1.23716634817820021358 * 0.01 + z * p2;
                p2 = 3.01581553508235416007 * 0.0001 + z * p2;
                p2 = 2.65806974686737550832 * 0.000001 + z * p2;
                p2 = 6.23974539184983293730 * 0.000000001 + z * p2;
                q2 = 1;
                q2 = 6.02427039364742014255 + z * q2;
                q2 = 3.67983563856160859403 + z * q2;
                q2 = 1.37702099489081330271 + z * q2;
                q2 = 2.16236993594496635890 * 0.1 + z * q2;
                q2 = 1.34204006088543189037 * 0.01 + z * q2;
                q2 = 3.28014464682127739104 * 0.0001 + z * q2;
                q2 = 2.89247864745380683936 * 0.000001 + z * q2;
                q2 = 6.79019408009981274425 * 0.000000001 + z * q2;
                x1 = z * p2 / q2;
            }
            x = x0 - x1;
            if (code != 0)
            {
                x = -x;
            }
            result = x;
            return result;
        }

        internal static int ablasblocksize(double[,] a)
        {
            int result = 0;

            result = 32;
            return result;
        }

        internal static void ablassplitlength(double[,] a,
   int n,
   ref int n1,
   ref int n2)
        {
            n1 = 0;
            n2 = 0;

            if (n > ablasblocksize(a))
            {
                ablasinternalsplitlength(n, ablasblocksize(a), ref n1, ref n2);
            }
            else
            {
                ablasinternalsplitlength(n, ablasmicroblocksize(), ref n1, ref n2);
            }
        }

        internal static void ablasinternalsplitlength(int n,
          int nb,
          ref int n1,
          ref int n2)
        {
            int r = 0;

            n1 = 0;
            n2 = 0;

            if (n <= nb)
            {

                //
                // Block size, no further splitting
                //
                n1 = n;
                n2 = 0;
            }
            else
            {

                //
                // Greater than block size
                //
                if (n % nb != 0)
                {

                    //
                    // Split remainder
                    //
                    n2 = n % nb;
                    n1 = n - n2;
                }
                else
                {

                    //
                    // Split on block boundaries
                    //
                    n2 = n / 2;
                    n1 = n - n2;
                    if (n1 % nb == 0)
                    {
                        return;
                    }
                    r = nb - n1 % nb;
                    n1 = n1 + r;
                    n2 = n2 - r;
                }
            }
        }

        internal static int ablasmicroblocksize()
        {
            int result = 0;

            result = 8;
            return result;
        }


        internal static void rmatrixgemm(int m,
   int n,
   int k,
   double alpha,
   double[,] a,
   int ia,
   int ja,
   int optypea,
   double[,] b,
   int ib,
   int jb,
   int optypeb,
   double beta,
   ref double[,] c,
   int ic,
   int jc)
        {
            int s1 = 0;
            int s2 = 0;
            int bs = 0;

            bs = ablasblocksize(a);
            if ((m <= bs & n <= bs) & k <= bs)
            {
                rmatrixgemmk(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                return;
            }
            if (m >= n & m >= k)
            {

                //
                // A*B = (A1 A2)^T*B
                //
                ablassplitlength(a, m, ref s1, ref s2);
                if (optypea == 0)
                {
                    rmatrixgemm(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(s2, n, k, alpha, a, ia + s1, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic + s1, jc);
                }
                else
                {
                    rmatrixgemm(s1, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(s2, n, k, alpha, a, ia, ja + s1, optypea, b, ib, jb, optypeb, beta, ref c, ic + s1, jc);
                }
                return;
            }
            if (n >= m & n >= k)
            {

                //
                // A*B = A*(B1 B2)
                //
                ablassplitlength(a, n, ref s1, ref s2);
                if (optypeb == 0)
                {
                    rmatrixgemm(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(m, s2, k, alpha, a, ia, ja, optypea, b, ib, jb + s1, optypeb, beta, ref c, ic, jc + s1);
                }
                else
                {
                    rmatrixgemm(m, s1, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(m, s2, k, alpha, a, ia, ja, optypea, b, ib + s1, jb, optypeb, beta, ref c, ic, jc + s1);
                }
                return;
            }
            if (k >= m & k >= n)
            {

                //
                // A*B = (A1 A2)*(B1 B2)^T
                //
                ablassplitlength(a, k, ref s1, ref s2);
                if (optypea == 0 & optypeb == 0)
                {
                    rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(m, n, s2, alpha, a, ia, ja + s1, optypea, b, ib + s1, jb, optypeb, 1.0, ref c, ic, jc);
                }
                if (optypea == 0 & optypeb != 0)
                {
                    rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(m, n, s2, alpha, a, ia, ja + s1, optypea, b, ib, jb + s1, optypeb, 1.0, ref c, ic, jc);
                }
                if (optypea != 0 & optypeb == 0)
                {
                    rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(m, n, s2, alpha, a, ia + s1, ja, optypea, b, ib + s1, jb, optypeb, 1.0, ref c, ic, jc);
                }
                if (optypea != 0 & optypeb != 0)
                {
                    rmatrixgemm(m, n, s1, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc);
                    rmatrixgemm(m, n, s2, alpha, a, ia + s1, ja, optypea, b, ib, jb + s1, optypeb, 1.0, ref c, ic, jc);
                }
                return;
            }
        }

        internal static void rmatrixgemmk(int m,
            int n,
            int k,
            double alpha,
            double[,] a,
            int ia,
            int ja,
            int optypea,
            double[,] b,
            int ib,
            int jb,
            int optypeb,
            double beta,
            ref double[,] c,
            int ic,
            int jc)
        {
            int i = 0;
            int j = 0;
            double v = 0;
            int i_ = 0;
            int i1_ = 0;


            //
            // if matrix size is zero
            //
            if (m * n == 0)
            {
                return;
            }

            //
            // Try optimized code
            //
            //if (rmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, ref c, ic, jc))
            //{
            //    return;
            //}

            //
            // if K=0, then C=Beta*C
            //
            if (k == 0)
            {
                if ((double)(beta) != (double)(1))
                {
                    if ((double)(beta) != (double)(0))
                    {
                        for (i = 0; i <= m - 1; i++)
                        {
                            for (j = 0; j <= n - 1; j++)
                            {
                                c[ic + i, jc + j] = beta * c[ic + i, jc + j];
                            }
                        }
                    }
                    else
                    {
                        for (i = 0; i <= m - 1; i++)
                        {
                            for (j = 0; j <= n - 1; j++)
                            {
                                c[ic + i, jc + j] = 0;
                            }
                        }
                    }
                }
                return;
            }

            //
            // General case
            //
            if (optypea == 0 & optypeb != 0)
            {

                //
                // A*B'
                //
                for (i = 0; i <= m - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        if (k == 0 | (double)(alpha) == (double)(0))
                        {
                            v = 0;
                        }
                        else
                        {
                            i1_ = (jb) - (ja);
                            v = 0.0;
                            for (i_ = ja; i_ <= ja + k - 1; i_++)
                            {
                                v += a[ia + i, i_] * b[ib + j, i_ + i1_];
                            }
                        }
                        if ((double)(beta) == (double)(0))
                        {
                            c[ic + i, jc + j] = alpha * v;
                        }
                        else
                        {
                            c[ic + i, jc + j] = beta * c[ic + i, jc + j] + alpha * v;
                        }
                    }
                }
                return;
            }
            if (optypea == 0 & optypeb == 0)
            {

                //
                // A*B
                //
                for (i = 0; i <= m - 1; i++)
                {
                    if ((double)(beta) != (double)(0))
                    {
                        for (i_ = jc; i_ <= jc + n - 1; i_++)
                        {
                            c[ic + i, i_] = beta * c[ic + i, i_];
                        }
                    }
                    else
                    {
                        for (j = 0; j <= n - 1; j++)
                        {
                            c[ic + i, jc + j] = 0;
                        }
                    }
                    if ((double)(alpha) != (double)(0))
                    {
                        for (j = 0; j <= k - 1; j++)
                        {
                            v = alpha * a[ia + i, ja + j];
                            i1_ = (jb) - (jc);
                            for (i_ = jc; i_ <= jc + n - 1; i_++)
                            {
                                c[ic + i, i_] = c[ic + i, i_] + v * b[ib + j, i_ + i1_];
                            }
                        }
                    }
                }
                return;
            }
            if (optypea != 0 & optypeb != 0)
            {

                //
                // A'*B'
                //
                for (i = 0; i <= m - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        if ((double)(alpha) == (double)(0))
                        {
                            v = 0;
                        }
                        else
                        {
                            i1_ = (jb) - (ia);
                            v = 0.0;
                            for (i_ = ia; i_ <= ia + k - 1; i_++)
                            {
                                v += a[i_, ja + i] * b[ib + j, i_ + i1_];
                            }
                        }
                        if ((double)(beta) == (double)(0))
                        {
                            c[ic + i, jc + j] = alpha * v;
                        }
                        else
                        {
                            c[ic + i, jc + j] = beta * c[ic + i, jc + j] + alpha * v;
                        }
                    }
                }
                return;
            }
            if (optypea != 0 & optypeb == 0)
            {

                //
                // A'*B
                //
                if ((double)(beta) == (double)(0))
                {
                    for (i = 0; i <= m - 1; i++)
                    {
                        for (j = 0; j <= n - 1; j++)
                        {
                            c[ic + i, jc + j] = 0;
                        }
                    }
                }
                else
                {
                    for (i = 0; i <= m - 1; i++)
                    {
                        for (i_ = jc; i_ <= jc + n - 1; i_++)
                        {
                            c[ic + i, i_] = beta * c[ic + i, i_];
                        }
                    }
                }
                if ((double)(alpha) != (double)(0))
                {
                    for (j = 0; j <= k - 1; j++)
                    {
                        for (i = 0; i <= m - 1; i++)
                        {
                            v = alpha * a[ia + j, ja + i];
                            i1_ = (jb) - (jc);
                            for (i_ = jc; i_ <= jc + n - 1; i_++)
                            {
                                c[ic + i, i_] = c[ic + i, i_] + v * b[ib + j, i_ + i1_];
                            }
                        }
                    }
                }
                return;
            }
        }


        public static double SpLine(List<KeyValuePair<double, double>> knownSamples, double z)
        {
            int np = knownSamples.Count;
            try
            {
                if (np > 1)
                {
                    double[] a = new double[np];
                    double x1;
                    double x2;
                    double y;
                    double[] h = new double[np];
                    for (int i = 1; i <= np - 1; i++)
                    {
                        h[i] = knownSamples[i].Key - knownSamples[i - 1].Key;
                    }
                    if (np > 2)
                    {
                        double[] sub = new double[np - 1];
                        double[] diag = new double[np - 1];
                        double[] sup = new double[np - 1];
                        for (int i = 1; i <= np - 2; i++)
                        {
                            diag[i] = (h[i] + h[i + 1]) / 3;
                            sup[i] = h[i + 1] / 6;
                            sub[i] = h[i] / 6;
                            a[i] = (knownSamples[i + 1].Value - knownSamples[i].Value) / h[i + 1] -
                                   (knownSamples[i].Value - knownSamples[i - 1].Value) / h[i];
                        }
                        solveTridiag(sub, diag, sup, ref a, np - 2);
                    }

                    int gap = 0;
                    //double previous = 0;
                    double previous = double.MinValue;
                    // At the end of this iteration, "gap" will contain the index of the interval
                    // between two known values, which contains the unknown z, and "previous" will
                    // contain the biggest z value among the known samples, left of the unknown z
                    for (int i = 0; i < knownSamples.Count; i++)
                    {
                        if (knownSamples[i].Key < z && knownSamples[i].Key > previous)
                        {
                            previous = knownSamples[i].Key;
                            gap = i + 1;
                        }
                    }
                    x1 = z - previous;
                    x2 = h[gap] - x1;
                    y = ((-a[gap - 1] / 6 * (x2 + h[gap]) * x1 + knownSamples[gap - 1].Value) * x2 +
                        (-a[gap] / 6 * (x1 + h[gap]) * x2 + knownSamples[gap].Value) * x1) / h[gap];
                    return y;
                }
            }
            catch (Exception)
            {

            }
            return 0;
        }

        static void solveTridiag(double[] sub, double[] diag, double[] sup, ref double[] b, int n)
        {

            int i;
            for (i = 2; i <= n; i++)
            {
                sub[i] = sub[i] / diag[i - 1];
                diag[i] = diag[i] - sub[i] * sup[i - 1];
                b[i] = b[i] - sub[i] * b[i - 1];
            }
            b[n] = b[n] / diag[n];
            for (i = n - 1; i >= 1; i--)
            {
                b[i] = (b[i] - sup[i] * b[i + 1]) / diag[i];
            }
        }

        public static bool spdmatrixcholesky(ref double[,] a,
            int n,
            bool isupper)
        {
            bool result = new bool();
            double[] tmp = new double[0];

            if (n < 1)
            {
                result = false;
                return result;
            }
            result = spdmatrixcholeskyrec(ref a, 0, n, isupper, ref tmp);
            return result;
        }

        public static bool spdmatrixcholeskyrec(ref double[,] a,
            int offs,
            int n,
            bool isupper,
            ref double[] tmp)
        {
            bool result = new bool();
            int n1 = 0;
            int n2 = 0;


            //
            // check N
            //
            if (n < 1)
            {
                result = false;
                return result;
            }

            //
            // Prepare buffer
            //
            if (tmp.Length < 2 * n)
            {
                tmp = new double[2 * n];
            }

            //
            // special cases
            //
            if (n == 1)
            {
                if ((double)(a[offs, offs]) > (double)(0))
                {
                    a[offs, offs] = Math.Sqrt(a[offs, offs]);
                    result = true;
                }
                else
                {
                    result = false;
                }
                return result;
            }
            if (n <= ablasblocksize(a))
            {
                result = spdmatrixcholesky2(ref a, offs, n, isupper, ref tmp);
                return result;
            }

            //
            // general case: split task in cache-oblivious manner
            //
            result = true;
            ablassplitlength(a, n, ref n1, ref n2);
            result = spdmatrixcholeskyrec(ref a, offs, n1, isupper, ref tmp);
            if (!result)
            {
                return result;
            }
            if (n2 > 0)
            {
                if (isupper)
                {
                    rmatrixlefttrsm(n1, n2, a, offs, offs, isupper, false, 1, ref a, offs, offs + n1);
                    rmatrixsyrk(n2, n1, -1.0, a, offs, offs + n1, 1, 1.0, ref a, offs + n1, offs + n1, isupper);
                }
                else
                {
                    rmatrixrighttrsm(n2, n1, a, offs, offs, isupper, false, 1, ref a, offs + n1, offs);
                    rmatrixsyrk(n2, n1, -1.0, a, offs + n1, offs, 0, 1.0, ref a, offs + n1, offs + n1, isupper);
                }
                result = spdmatrixcholeskyrec(ref a, offs + n1, n2, isupper, ref tmp);
                if (!result)
                {
                    return result;
                }
            }
            return result;
        }

    



        private static bool spdmatrixcholesky2(ref double[,] aaa,
            int offs,
            int n,
            bool isupper,
            ref double[] tmp)
        {
            bool result = new bool();
            int i = 0;
            int j = 0;
            double ajj = 0;
            double v = 0;
            double r = 0;
            int i_ = 0;
            int i1_ = 0;

            result = true;
            if (n < 0)
            {
                result = false;
                return result;
            }

            //
            // Quick return if possible
            //
            if (n == 0)
            {
                return result;
            }
            if (isupper)
            {

                //
                // Compute the Cholesky factorization A = U'*U.
                //
                for (j = 0; j <= n - 1; j++)
                {

                    //
                    // Compute U(J,J) and test for non-positive-definiteness.
                    //
                    v = 0.0;
                    for (i_ = offs; i_ <= offs + j - 1; i_++)
                    {
                        v += aaa[i_, offs + j] * aaa[i_, offs + j];
                    }
                    ajj = aaa[offs + j, offs + j] - v;
                    if ((double)(ajj) <= (double)(0))
                    {
                        aaa[offs + j, offs + j] = ajj;
                        result = false;
                        return result;
                    }
                    ajj = Math.Sqrt(ajj);
                    aaa[offs + j, offs + j] = ajj;

                    //
                    // Compute elements J+1:N-1 of row J.
                    //
                    if (j < n - 1)
                    {
                        if (j > 0)
                        {
                            i1_ = (offs) - (0);
                            for (i_ = 0; i_ <= j - 1; i_++)
                            {
                                tmp[i_] = -aaa[i_ + i1_, offs + j];
                            }
                            rmatrixmv(n - j - 1, j, aaa, offs, offs + j + 1, 1, tmp, 0, ref tmp, n);
                            i1_ = (n) - (offs + j + 1);
                            for (i_ = offs + j + 1; i_ <= offs + n - 1; i_++)
                            {
                                aaa[offs + j, i_] = aaa[offs + j, i_] + tmp[i_ + i1_];
                            }
                        }
                        r = 1 / ajj;
                        for (i_ = offs + j + 1; i_ <= offs + n - 1; i_++)
                        {
                            aaa[offs + j, i_] = r * aaa[offs + j, i_];
                        }
                    }
                }
            }
            else
            {

                //
                // Compute the Cholesky factorization A = L*L'.
                //
                for (j = 0; j <= n - 1; j++)
                {

                    //
                    // Compute L(J+1,J+1) and test for non-positive-definiteness.
                    //
                    v = 0.0;
                    for (i_ = offs; i_ <= offs + j - 1; i_++)
                    {
                        v += aaa[offs + j, i_] * aaa[offs + j, i_];
                    }
                    ajj = aaa[offs + j, offs + j] - v;
                    if ((double)(ajj) <= (double)(0))
                    {
                        aaa[offs + j, offs + j] = ajj;
                        result = false;
                        return result;
                    }
                    ajj = Math.Sqrt(ajj);
                    aaa[offs + j, offs + j] = ajj;

                    //
                    // Compute elements J+1:N of column J.
                    //
                    if (j < n - 1)
                    {
                        if (j > 0)
                        {
                            i1_ = (offs) - (0);
                            for (i_ = 0; i_ <= j - 1; i_++)
                            {
                                tmp[i_] = aaa[offs + j, i_ + i1_];
                            }
                            rmatrixmv(n - j - 1, j, aaa, offs + j + 1, offs, 0, tmp, 0, ref tmp, n);
                            for (i = 0; i <= n - j - 2; i++)
                            {
                                aaa[offs + j + 1 + i, offs + j] = (aaa[offs + j + 1 + i, offs + j] - tmp[n + i]) / ajj;
                            }
                        }
                        else
                        {
                            for (i = 0; i <= n - j - 2; i++)
                            {
                                aaa[offs + j + 1 + i, offs + j] = aaa[offs + j + 1 + i, offs + j] / ajj;
                            }
                        }
                    }
                }
            }
            return result;
        }

        public static void rmatrixmv(int m,
           int n,
           double[,] a,
           int ia,
           int ja,
           int opa,
           double[] x,
           int ix,
           ref double[] y,
           int iy)
        {
            int i = 0;
            double v = 0;
            int i_ = 0;
            int i1_ = 0;

            if (m == 0)
            {
                return;
            }
            if (n == 0)
            {
                for (i = 0; i <= m - 1; i++)
                {
                    y[iy + i] = 0;
                }
                return;
            }
            if (rmatrixmvf(m, n, a, ia, ja, opa, x, ix, ref y, iy))
            {
                return;
            }
            if (opa == 0)
            {

                //
                // y = A*x
                //
                for (i = 0; i <= m - 1; i++)
                {
                    i1_ = (ix) - (ja);
                    v = 0.0;
                    for (i_ = ja; i_ <= ja + n - 1; i_++)
                    {
                        v += a[ia + i, i_] * x[i_ + i1_];
                    }
                    y[iy + i] = v;
                }
                return;
            }
            if (opa == 1)
            {

                //
                // y = A^T*x
                //
                for (i = 0; i <= m - 1; i++)
                {
                    y[iy + i] = 0;
                }
                for (i = 0; i <= n - 1; i++)
                {
                    v = x[ix + i];
                    i1_ = (ja) - (iy);
                    for (i_ = iy; i_ <= iy + m - 1; i_++)
                    {
                        y[i_] = y[i_] + v * a[ia + i, i_ + i1_];
                    }
                }
                return;
            }
        }

        public static void rmatrixsyrk(int n,
             int k,
             double alpha,
             double[,] a,
             int ia,
             int ja,
             int optypea,
             double beta,
             ref double[,] c,
             int ic,
             int jc,
             bool isupper)
        {
            int s1 = 0;
            int s2 = 0;
            int bs = 0;

            bs = ablasblocksize(a);
            if (n <= bs & k <= bs)
            {
                rmatrixsyrk2(n, k, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);
                return;
            }
            if (k >= n)
            {

                //
                // Split K
                //
                ablassplitlength(a, k, ref s1, ref s2);
                if (optypea == 0)
                {
                    rmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);
                    rmatrixsyrk(n, s2, alpha, a, ia, ja + s1, optypea, 1.0, ref c, ic, jc, isupper);
                }
                else
                {
                    rmatrixsyrk(n, s1, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);
                    rmatrixsyrk(n, s2, alpha, a, ia + s1, ja, optypea, 1.0, ref c, ic, jc, isupper);
                }
            }
            else
            {

                //
                // Split N
                //
                ablassplitlength(a, n, ref s1, ref s2);
                if (optypea == 0 & isupper)
                {
                    rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);
                    rmatrixgemm(s1, s2, k, alpha, a, ia, ja, 0, a, ia + s1, ja, 1, beta, ref c, ic, jc + s1);
                    rmatrixsyrk(s2, k, alpha, a, ia + s1, ja, optypea, beta, ref c, ic + s1, jc + s1, isupper);
                    return;
                }
                if (optypea == 0 & !isupper)
                {
                    rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);
                    rmatrixgemm(s2, s1, k, alpha, a, ia + s1, ja, 0, a, ia, ja, 1, beta, ref c, ic + s1, jc);
                    rmatrixsyrk(s2, k, alpha, a, ia + s1, ja, optypea, beta, ref c, ic + s1, jc + s1, isupper);
                    return;
                }
                if (optypea != 0 & isupper)
                {
                    rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);
                    rmatrixgemm(s1, s2, k, alpha, a, ia, ja, 1, a, ia, ja + s1, 0, beta, ref c, ic, jc + s1);
                    rmatrixsyrk(s2, k, alpha, a, ia, ja + s1, optypea, beta, ref c, ic + s1, jc + s1, isupper);
                    return;
                }
                if (optypea != 0 & !isupper)
                {
                    rmatrixsyrk(s1, k, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper);
                    rmatrixgemm(s2, s1, k, alpha, a, ia, ja + s1, 1, a, ia, ja, 0, beta, ref c, ic + s1, jc);
                    rmatrixsyrk(s2, k, alpha, a, ia, ja + s1, optypea, beta, ref c, ic + s1, jc + s1, isupper);
                    return;
                }
            }
        }

        private static void rmatrixsyrk2(int n,
           int k,
           double alpha,
           double[,] a,
           int ia,
           int ja,
           int optypea,
           double beta,
           ref double[,] c,
           int ic,
           int jc,
           bool isupper)
        {
            int i = 0;
            int j = 0;
            int j1 = 0;
            int j2 = 0;
            double v = 0;
            int i_ = 0;
            int i1_ = 0;


            //
            // Fast exit (nothing to be done)
            //
            if (((double)(alpha) == (double)(0) | k == 0) & (double)(beta) == (double)(1))
            {
                return;
            }

            //
            // Try to call fast SYRK
            //
            //if (rmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, ref c, ic, jc, isupper))
            //{
            //    return;
            //}

            //
            // SYRK
            //
            if (optypea == 0)
            {

                //
                // C=alpha*A*A^H+beta*C
                //
                for (i = 0; i <= n - 1; i++)
                {
                    if (isupper)
                    {
                        j1 = i;
                        j2 = n - 1;
                    }
                    else
                    {
                        j1 = 0;
                        j2 = i;
                    }
                    for (j = j1; j <= j2; j++)
                    {
                        if ((double)(alpha) != (double)(0) & k > 0)
                        {
                            v = 0.0;
                            for (i_ = ja; i_ <= ja + k - 1; i_++)
                            {
                                v += a[ia + i, i_] * a[ia + j, i_];
                            }
                        }
                        else
                        {
                            v = 0;
                        }
                        if ((double)(beta) == (double)(0))
                        {
                            c[ic + i, jc + j] = alpha * v;
                        }
                        else
                        {
                            c[ic + i, jc + j] = beta * c[ic + i, jc + j] + alpha * v;
                        }
                    }
                }
                return;
            }
            else
            {

                //
                // C=alpha*A^H*A+beta*C
                //
                for (i = 0; i <= n - 1; i++)
                {
                    if (isupper)
                    {
                        j1 = i;
                        j2 = n - 1;
                    }
                    else
                    {
                        j1 = 0;
                        j2 = i;
                    }
                    if ((double)(beta) == (double)(0))
                    {
                        for (j = j1; j <= j2; j++)
                        {
                            c[ic + i, jc + j] = 0;
                        }
                    }
                    else
                    {
                        for (i_ = jc + j1; i_ <= jc + j2; i_++)
                        {
                            c[ic + i, i_] = beta * c[ic + i, i_];
                        }
                    }
                }
                for (i = 0; i <= k - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        if (isupper)
                        {
                            j1 = j;
                            j2 = n - 1;
                        }
                        else
                        {
                            j1 = 0;
                            j2 = j;
                        }
                        v = alpha * a[ia + i, ja + j];
                        i1_ = (ja + j1) - (jc + j1);
                        for (i_ = jc + j1; i_ <= jc + j2; i_++)
                        {
                            c[ic + j, i_] = c[ic + j, i_] + v * a[ia + i, i_ + i1_];
                        }
                    }
                }
                return;
            }
        }

        public static void rmatrixlefttrsm(int m,
           int n,
           double[,] a,
           int i1,
           int j1,
           bool isupper,
           bool isunit,
           int optype,
           ref double[,] x,
           int i2,
           int j2)
        {
            int s1 = 0;
            int s2 = 0;
            int bs = 0;

            bs = ablasblocksize(a);
            if (m <= bs & n <= bs)
            {
                rmatrixlefttrsm2(m, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                return;
            }
            if (n >= m)
            {

                //
                // Split X: op(A)^-1*X = op(A)^-1*(X1 X2)
                //
                ablassplitlength(x, n, ref s1, ref s2);
                rmatrixlefttrsm(m, s1, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                rmatrixlefttrsm(m, s2, a, i1, j1, isupper, isunit, optype, ref x, i2, j2 + s1);
            }
            else
            {

                //
                // Split A
                //
                ablassplitlength(a, m, ref s1, ref s2);
                if (isupper & optype == 0)
                {

                    //
                    //           (A1  A12)-1  ( X1 )
                    // A^-1*X* = (       )   *(    )
                    //           (     A2)    ( X2 )
                    //
                    rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2 + s1, j2);
                    rmatrixgemm(s1, n, s2, -1.0, a, i1, j1 + s1, 0, x, i2 + s1, j2, 0, 1.0, ref x, i2, j2);
                    rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    return;
                }
                if (isupper & optype != 0)
                {

                    //
                    //          (A1'     )-1 ( X1 )
                    // A^-1*X = (        )  *(    )
                    //          (A12' A2')   ( X2 )
                    //
                    rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    rmatrixgemm(s2, n, s1, -1.0, a, i1, j1 + s1, optype, x, i2, j2, 0, 1.0, ref x, i2 + s1, j2);
                    rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2 + s1, j2);
                    return;
                }
                if (!isupper & optype == 0)
                {

                    //
                    //          (A1     )-1 ( X1 )
                    // A^-1*X = (       )  *(    )
                    //          (A21  A2)   ( X2 )
                    //
                    rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    rmatrixgemm(s2, n, s1, -1.0, a, i1 + s1, j1, 0, x, i2, j2, 0, 1.0, ref x, i2 + s1, j2);
                    rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2 + s1, j2);
                    return;
                }
                if (!isupper & optype != 0)
                {

                    //
                    //          (A1' A21')-1 ( X1 )
                    // A^-1*X = (        )  *(    )
                    //          (     A2')   ( X2 )
                    //
                    rmatrixlefttrsm(s2, n, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2 + s1, j2);
                    rmatrixgemm(s1, n, s2, -1.0, a, i1 + s1, j1, optype, x, i2 + s1, j2, 0, 1.0, ref x, i2, j2);
                    rmatrixlefttrsm(s1, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    return;
                }
            }
        }

        private static void rmatrixlefttrsm2(int m,
           int n,
           double[,] a,
           int i1,
           int j1,
           bool isupper,
           bool isunit,
           int optype,
           ref double[,] x,
           int i2,
           int j2)
        {
            int i = 0;
            int j = 0;
            double vr = 0;
            double vd = 0;
            int i_ = 0;


            //
            // Special case
            //
            if (n * m == 0)
            {
                return;
            }

            //
            // Try fast code
            //
            //if (rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2))
            //{
            //    return;
            //}

            //
            // General case
            //
            if (isupper)
            {

                //
                // Upper triangular matrix
                //
                if (optype == 0)
                {

                    //
                    // A^(-1)*X
                    //
                    for (i = m - 1; i >= 0; i--)
                    {
                        for (j = i + 1; j <= m - 1; j++)
                        {
                            vr = a[i1 + i, j1 + j];
                            for (i_ = j2; i_ <= j2 + n - 1; i_++)
                            {
                                x[i2 + i, i_] = x[i2 + i, i_] - vr * x[i2 + j, i_];
                            }
                        }
                        if (!isunit)
                        {
                            vd = 1 / a[i1 + i, j1 + i];
                            for (i_ = j2; i_ <= j2 + n - 1; i_++)
                            {
                                x[i2 + i, i_] = vd * x[i2 + i, i_];
                            }
                        }
                    }
                    return;
                }
                if (optype == 1)
                {

                    //
                    // A^(-T)*X
                    //
                    for (i = 0; i <= m - 1; i++)
                    {
                        if (isunit)
                        {
                            vd = 1;
                        }
                        else
                        {
                            vd = 1 / a[i1 + i, j1 + i];
                        }
                        for (i_ = j2; i_ <= j2 + n - 1; i_++)
                        {
                            x[i2 + i, i_] = vd * x[i2 + i, i_];
                        }
                        for (j = i + 1; j <= m - 1; j++)
                        {
                            vr = a[i1 + i, j1 + j];
                            for (i_ = j2; i_ <= j2 + n - 1; i_++)
                            {
                                x[i2 + j, i_] = x[i2 + j, i_] - vr * x[i2 + i, i_];
                            }
                        }
                    }
                    return;
                }
            }
            else
            {

                //
                // Lower triangular matrix
                //
                if (optype == 0)
                {

                    //
                    // A^(-1)*X
                    //
                    for (i = 0; i <= m - 1; i++)
                    {
                        for (j = 0; j <= i - 1; j++)
                        {
                            vr = a[i1 + i, j1 + j];
                            for (i_ = j2; i_ <= j2 + n - 1; i_++)
                            {
                                x[i2 + i, i_] = x[i2 + i, i_] - vr * x[i2 + j, i_];
                            }
                        }
                        if (isunit)
                        {
                            vd = 1;
                        }
                        else
                        {
                            vd = 1 / a[i1 + j, j1 + j];
                        }
                        for (i_ = j2; i_ <= j2 + n - 1; i_++)
                        {
                            x[i2 + i, i_] = vd * x[i2 + i, i_];
                        }
                    }
                    return;
                }
                if (optype == 1)
                {

                    //
                    // A^(-T)*X
                    //
                    for (i = m - 1; i >= 0; i--)
                    {
                        if (isunit)
                        {
                            vd = 1;
                        }
                        else
                        {
                            vd = 1 / a[i1 + i, j1 + i];
                        }
                        for (i_ = j2; i_ <= j2 + n - 1; i_++)
                        {
                            x[i2 + i, i_] = vd * x[i2 + i, i_];
                        }
                        for (j = i - 1; j >= 0; j--)
                        {
                            vr = a[i1 + i, j1 + j];
                            for (i_ = j2; i_ <= j2 + n - 1; i_++)
                            {
                                x[i2 + j, i_] = x[i2 + j, i_] - vr * x[i2 + i, i_];
                            }
                        }
                    }
                    return;
                }
            }
        }

        public static void rmatrixrighttrsm(int m,
           int n,
           double[,] a,
           int i1,
           int j1,
           bool isupper,
           bool isunit,
           int optype,
           ref double[,] x,
           int i2,
           int j2)
        {
            int s1 = 0;
            int s2 = 0;
            int bs = 0;

            bs = ablasblocksize(a);
            if (m <= bs & n <= bs)
            {
                rmatrixrighttrsm2(m, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                return;
            }
            if (m >= n)
            {

                //
                // Split X: X*A = (X1 X2)^T*A
                //
                ablassplitlength(a, m, ref s1, ref s2);
                rmatrixrighttrsm(s1, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                rmatrixrighttrsm(s2, n, a, i1, j1, isupper, isunit, optype, ref x, i2 + s1, j2);
            }
            else
            {

                //
                // Split A:
                //               (A1  A12)
                // X*op(A) = X*op(       )
                //               (     A2)
                //
                // Different variants depending on
                // IsUpper/OpType combinations
                //
                ablassplitlength(a, n, ref s1, ref s2);
                if (isupper & optype == 0)
                {

                    //
                    //                  (A1  A12)-1
                    // X*A^-1 = (X1 X2)*(       )
                    //                  (     A2)
                    //
                    rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    rmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1, j1 + s1, 0, 1.0, ref x, i2, j2 + s1);
                    rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2, j2 + s1);
                    return;
                }
                if (isupper & optype != 0)
                {

                    //
                    //                  (A1'     )-1
                    // X*A^-1 = (X1 X2)*(        )
                    //                  (A12' A2')
                    //
                    rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2, j2 + s1);
                    rmatrixgemm(m, s1, s2, -1.0, x, i2, j2 + s1, 0, a, i1, j1 + s1, optype, 1.0, ref x, i2, j2);
                    rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    return;
                }
                if (!isupper & optype == 0)
                {

                    //
                    //                  (A1     )-1
                    // X*A^-1 = (X1 X2)*(       )
                    //                  (A21  A2)
                    //
                    rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2, j2 + s1);
                    rmatrixgemm(m, s1, s2, -1.0, x, i2, j2 + s1, 0, a, i1 + s1, j1, 0, 1.0, ref x, i2, j2);
                    rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    return;
                }
                if (!isupper & optype != 0)
                {

                    //
                    //                  (A1' A21')-1
                    // X*A^-1 = (X1 X2)*(        )
                    //                  (     A2')
                    //
                    rmatrixrighttrsm(m, s1, a, i1, j1, isupper, isunit, optype, ref x, i2, j2);
                    rmatrixgemm(m, s2, s1, -1.0, x, i2, j2, 0, a, i1 + s1, j1, optype, 1.0, ref x, i2, j2 + s1);
                    rmatrixrighttrsm(m, s2, a, i1 + s1, j1 + s1, isupper, isunit, optype, ref x, i2, j2 + s1);
                    return;
                }
            }
        }

        private static void rmatrixrighttrsm2(int m,
            int n,
            double[,] a,
            int i1,
            int j1,
            bool isupper,
            bool isunit,
            int optype,
            ref double[,] x,
            int i2,
            int j2)
        {
            int i = 0;
            int j = 0;
            double vr = 0;
            double vd = 0;
            int i_ = 0;
            int i1_ = 0;


            //
            // Special case
            //
            if (n * m == 0)
            {
                return;
            }

            //
            // Try to use "fast" code
            //
            if (rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, ref x, i2, j2))
            {
                return;
            }

            //
            // General case
            //
            if (isupper)
            {

                //
                // Upper triangular matrix
                //
                if (optype == 0)
                {

                    //
                    // X*A^(-1)
                    //
                    for (i = 0; i <= m - 1; i++)
                    {
                        for (j = 0; j <= n - 1; j++)
                        {
                            if (isunit)
                            {
                                vd = 1;
                            }
                            else
                            {
                                vd = a[i1 + j, j1 + j];
                            }
                            x[i2 + i, j2 + j] = x[i2 + i, j2 + j] / vd;
                            if (j < n - 1)
                            {
                                vr = x[i2 + i, j2 + j];
                                i1_ = (j1 + j + 1) - (j2 + j + 1);
                                for (i_ = j2 + j + 1; i_ <= j2 + n - 1; i_++)
                                {
                                    x[i2 + i, i_] = x[i2 + i, i_] - vr * a[i1 + j, i_ + i1_];
                                }
                            }
                        }
                    }
                    return;
                }
                if (optype == 1)
                {

                    //
                    // X*A^(-T)
                    //
                    for (i = 0; i <= m - 1; i++)
                    {
                        for (j = n - 1; j >= 0; j--)
                        {
                            vr = 0;
                            vd = 1;
                            if (j < n - 1)
                            {
                                i1_ = (j1 + j + 1) - (j2 + j + 1);
                                vr = 0.0;
                                for (i_ = j2 + j + 1; i_ <= j2 + n - 1; i_++)
                                {
                                    vr += x[i2 + i, i_] * a[i1 + j, i_ + i1_];
                                }
                            }
                            if (!isunit)
                            {
                                vd = a[i1 + j, j1 + j];
                            }
                            x[i2 + i, j2 + j] = (x[i2 + i, j2 + j] - vr) / vd;
                        }
                    }
                    return;
                }
            }
            else
            {

                //
                // Lower triangular matrix
                //
                if (optype == 0)
                {

                    //
                    // X*A^(-1)
                    //
                    for (i = 0; i <= m - 1; i++)
                    {
                        for (j = n - 1; j >= 0; j--)
                        {
                            if (isunit)
                            {
                                vd = 1;
                            }
                            else
                            {
                                vd = a[i1 + j, j1 + j];
                            }
                            x[i2 + i, j2 + j] = x[i2 + i, j2 + j] / vd;
                            if (j > 0)
                            {
                                vr = x[i2 + i, j2 + j];
                                i1_ = (j1) - (j2);
                                for (i_ = j2; i_ <= j2 + j - 1; i_++)
                                {
                                    x[i2 + i, i_] = x[i2 + i, i_] - vr * a[i1 + j, i_ + i1_];
                                }
                            }
                        }
                    }
                    return;
                }
                if (optype == 1)
                {

                    //
                    // X*A^(-T)
                    //
                    for (i = 0; i <= m - 1; i++)
                    {
                        for (j = 0; j <= n - 1; j++)
                        {
                            vr = 0;
                            vd = 1;
                            if (j > 0)
                            {
                                i1_ = (j1) - (j2);
                                vr = 0.0;
                                for (i_ = j2; i_ <= j2 + j - 1; i_++)
                                {
                                    vr += x[i2 + i, i_] * a[i1 + j, i_ + i1_];
                                }
                            }
                            if (!isunit)
                            {
                                vd = a[i1 + j, j1 + j];
                            }
                            x[i2 + i, j2 + j] = (x[i2 + i, j2 + j] - vr) / vd;
                        }
                    }
                    return;
                }
            }
        }

        public static bool rmatrixrighttrsmf(int m,
           int n,
           double[,] a,
           int i1,
           int j1,
           bool isupper,
           bool isunit,
           int optype,
           ref double[,] x,
           int i2,
           int j2)
        {
            bool result = new bool();

            result = false;
            return result;
        }

        public static bool rmatrixmvf(int m,
           int n,
           double[,] a,
           int ia,
           int ja,
           int opa,
           double[] x,
           int ix,
           ref double[] y,
           int iy)
        {
            bool result = new bool();

            result = false;
            return result;
        }
    }
}

public class SimpleRNG
{
    private static uint m_w;
    private static uint m_z;

    static SimpleRNG()
    {
        // These values are not magical, just the default values Marsaglia used.
        // Any pair of unsigned integers should be fine.
        m_w = 521288629;
        m_z = 362436069;
    }

    // The random generator seed can be set three ways:
    // 1) specifying two non-zero unsigned integers
    // 2) specifying one non-zero unsigned integer and taking a default value for the second
    // 3) setting the seed from the system time

    public static void SetSeed(uint u, uint v)
    {
        if (u != 0) m_w = u;
        if (v != 0) m_z = v;
    }

    public static void SetSeed(uint u)
    {
        m_w = u;
    }

    public static void SetSeedFromSystemTime()
    {
        System.DateTime dt = System.DateTime.Now;
        long x = dt.ToFileTime();
        SetSeed((uint)(x >> 16), (uint)(x % 4294967296));
    }

    // Produce a uniform random sample from the open interval (0, 1).
    // The method will not return either end point.
    public static double GetUniform()
    {
        // 0 <= u < 2^32
        uint u = GetUint();
        // The magic number below is 1/(2^32 + 2).
        // The result is strictly between 0 and 1.
        return (u + 1.0) * 2.328306435454494e-10;
    }

    // This is the heart of the generator.
    // It uses George Marsaglia's MWC algorithm to produce an unsigned integer.
    // See http://www.bobwheeler.com/statistics/Password/MarsagliaPost.txt
    private static uint GetUint()
    {
        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
        return (m_z << 16) + m_w;
    }

    // Get normal (Gaussian) random sample with mean 0 and standard deviation 1
    public static double GetNormal()
    {
        // Use Box-Muller algorithm
        double u1 = GetUniform();
        double u2 = GetUniform();
        double r = Math.Sqrt(-2.0 * Math.Log(u1));
        double theta = 2.0 * Math.PI * u2;
        return r * Math.Sin(theta);
    }

    // Get normal (Gaussian) random sample with specified mean and standard deviation
    public static double GetNormal(double mean, double standardDeviation)
    {
        if (standardDeviation <= 0.0)
        {
            string msg = string.Format("Shape must be positive. Received {0}.", standardDeviation);
            throw new ArgumentOutOfRangeException(msg);
        }
        return mean + standardDeviation * GetNormal();
    }

    // Get exponential random sample with mean 1
    public static double GetExponential()
    {
        return -Math.Log(GetUniform());
    }

    // Get exponential random sample with specified mean
    public static double GetExponential(double mean)
    {
        if (mean <= 0.0)
        {
            string msg = string.Format("Mean must be positive. Received {0}.", mean);
            throw new ArgumentOutOfRangeException(msg);
        }
        return mean * GetExponential();
    }

    public static double GetGamma(double shape, double scale)
    {
        // Implementation based on "A Simple Method for Generating Gamma Variables"
        // by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
        // Vol 26, No 3, September 2000, pages 363-372.

        double d, c, x, xsquared, v, u;

        if (shape >= 1.0)
        {
            d = shape - 1.0 / 3.0;
            c = 1.0 / Math.Sqrt(9.0 * d);
            for (;;)
            {
                do
                {
                    x = GetNormal();
                    v = 1.0 + c * x;
                }
                while (v <= 0.0);
                v = v * v * v;
                u = GetUniform();
                xsquared = x * x;
                if (u < 1.0 - .0331 * xsquared * xsquared || Math.Log(u) < 0.5 * xsquared + d * (1.0 - v + Math.Log(v)))
                    return scale * d * v;
            }
        }
        else if (shape <= 0.0)
        {
            string msg = string.Format("Shape must be positive. Received {0}.", shape);
            throw new ArgumentOutOfRangeException(msg);
        }
        else
        {
            double g = GetGamma(shape + 1.0, 1.0);
            double w = GetUniform();
            return scale * g * Math.Pow(w, 1.0 / shape);
        }
    }



}


public class SimpleRNGNew
{
    private static uint m_w;
    private static uint m_z;

    public SimpleRNGNew()
    {
        // These values are not magical, just the default values Marsaglia used.
        // Any pair of unsigned integers should be fine.
        m_w = 521288629;
        m_z = 362436069;
    }

    // The random generator seed can be set three ways:
    // 1) specifying two non-zero unsigned integers
    // 2) specifying one non-zero unsigned integer and taking a default value for the second
    // 3) setting the seed from the system time

    public void SetSeed(uint u, uint v)
    {
        if (u != 0) m_w = u;
        if (v != 0) m_z = v;
    }

    public void SetSeed(uint u)
    {
        m_w = u;
    }

    public void SetSeedFromSystemTime()
    {
        System.DateTime dt = System.DateTime.Now;
        long x = dt.ToFileTime();
        SetSeed((uint)(x >> 16), (uint)(x % 4294967296));
    }

    // Produce a uniform random sample from the open interval (0, 1).
    // The method will not return either end point.
    public double GetUniform()
    {
        // 0 <= u < 2^32
        uint u = GetUint();
        // The magic number below is 1/(2^32 + 2).
        // The result is strictly between 0 and 1.
        return (u + 1.0) * 2.328306435454494e-10;
    }

    // This is the heart of the generator.
    // It uses George Marsaglia's MWC algorithm to produce an unsigned integer.
    // See http://www.bobwheeler.com/statistics/Password/MarsagliaPost.txt
    private uint GetUint()
    {
        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
        return (m_z << 16) + m_w;
    }

    // Get normal (Gaussian) random sample with mean 0 and standard deviation 1
    public double GetNormal()
    {
        // Use Box-Muller algorithm
        double u1 = GetUniform();
        double u2 = GetUniform();
        double r = Math.Sqrt(-2.0 * Math.Log(u1));
        double theta = 2.0 * Math.PI * u2;
        return r * Math.Sin(theta);
    }

    // Get normal (Gaussian) random sample with specified mean and standard deviation
    public double GetNormal(double mean, double standardDeviation)
    {
        if (standardDeviation <= 0.0)
        {
            string msg = string.Format("Shape must be positive. Received {0}.", standardDeviation);
            throw new ArgumentOutOfRangeException(msg);
        }
        return mean + standardDeviation * GetNormal();
    }

    // Get exponential random sample with mean 1
    public double GetExponential()
    {
        return -Math.Log(GetUniform());
    }

    // Get exponential random sample with specified mean
    public double GetExponential(double mean)
    {
        if (mean <= 0.0)
        {
            string msg = string.Format("Mean must be positive. Received {0}.", mean);
            throw new ArgumentOutOfRangeException(msg);
        }
        return mean * GetExponential();
    }

    public double GetGamma(double shape, double scale)
    {
        // Implementation based on "A Simple Method for Generating Gamma Variables"
        // by George Marsaglia and Wai Wan Tsang.  ACM Transactions on Mathematical Software
        // Vol 26, No 3, September 2000, pages 363-372.

        double d, c, x, xsquared, v, u;

        if (shape >= 1.0)
        {
            d = shape - 1.0 / 3.0;
            c = 1.0 / Math.Sqrt(9.0 * d);
            for (;;)
            {
                do
                {
                    x = GetNormal();
                    v = 1.0 + c * x;
                }
                while (v <= 0.0);
                v = v * v * v;
                u = GetUniform();
                xsquared = x * x;
                if (u < 1.0 - .0331 * xsquared * xsquared || Math.Log(u) < 0.5 * xsquared + d * (1.0 - v + Math.Log(v)))
                    return scale * d * v;
            }
        }
        else if (shape <= 0.0)
        {
            string msg = string.Format("Shape must be positive. Received {0}.", shape);
            throw new ArgumentOutOfRangeException(msg);
        }
        else
        {
            double g = GetGamma(shape + 1.0, 1.0);
            double w = GetUniform();
            return scale * g * Math.Pow(w, 1.0 / shape);
        }
    }
}


