using System;
using Redzen.Numerics.Distributions;
using Redzen.Numerics.Distributions.Double;
using Redzen.Random;

namespace FractionalBrownianMotion
{
    public class HoskingMethod
    {
        #region Instance Fields

        // Method parameters.
        readonly int _m;
        readonly double _h;
        readonly double _l;
        readonly double _scaling;

        // Working memory.
        readonly double[] _phi;
        readonly double[] _psi;
        readonly double[] _cov;

        readonly ISampler<double> _gaussian;

        #endregion

        #region Constructors

        /// <summary>
        /// 
        /// </summary>
        /// <param name="n">Determines the sample size N by N=2^(*n)</param>
        /// <param name="h">The Hurst parameter of the trace</param>
        /// <param name="l">The sample is generated on [0,L]</param>
        public HoskingMethod(int n, double h, double l)
        {
            _m = 1 << n;  // m = 2^n
            _h = h;
            _l = l;
            _scaling = Math.Pow(_l / _m, h);

            _phi = new double[_m];
            _psi = new double[_m];
            _cov = new double[_m];

            _gaussian = new ZigguratGaussianSampler();
        }

        public HoskingMethod(int n, double h, double l, ulong seed)
        {
            _m = 1 << n;  // m = 2^n
            _h = h;
            _l = l;
            _scaling = Math.Pow(_l / _m, h);

            _phi = new double[_m];
            _psi = new double[_m];
            _cov = new double[_m];

            _gaussian = new ZigguratGaussianSampler(0.0, 1.0, RandomDefaults.CreateRandomSource(seed));
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// Generates a set of fractional Brownian noise samples using the Hosking method.      
        /// </summary>
        /// <param name="output">An array to fill with samples.</param>
        public void Sample(double[] output)
        {           
            if(output.Length != _m) throw new ArgumentException("Invalid length");

            // Initialization.
            output[0] = _gaussian.Sample();
            double v = 1;
            _phi[0] = 0;
            for(int i=0; i < _m; i++)
            _cov[i] = Covariance(i, _h);

            // Simulation.
            for(int i=1; i < _m; i++) 
            {
                _phi[i-1] = _cov[i];
                for(int j=0; j < i-1; j++) 
                {
                    _psi[j] = _phi[j];
                    _phi[i-1] -= _psi[j] * _cov[i-j-1];
                }

                _phi[i-1] /= v;
                for(int j=0; j < i-1; j++) {
                    _phi[j] = _psi[j] - _phi[i-1] * _psi[i-j-2];
                }

                v *= (1 - _phi[i-1] * _phi[i-1]);
    
                output[i] = 0.0;
                for(int j=0; j < i; j++) {
                    output[i] += _phi[j] * output[i-j-1];
                }
                output[i] += Math.Sqrt(v) * _gaussian.Sample();
            }

            // Rescale to obtain a sample of size 2^(*n) on [0,L].
            for(int i=0; i < _m; i++) {
                output[i] = _scaling * (output[i]);
            }
        }

        #endregion

        #region Private Static Methods

        private static double Covariance(long i, double H) 
        {
            if (i == 0) return 1;
            else return (Math.Pow(i-1, 2.0*H) - 2*Math.Pow(i, 2.0*H) + Math.Pow(i+1, 2.0*H)) / 2.0;
        }

        #endregion
    }
}
