using System;
using Xunit;

namespace FractionalBrownianMotion.Tests
{
    public class HoskingMethodTests
    {
        [Fact]
        public void SampleNoise()
        {
            var hosking = new HoskingMethod(8, 0.7, 100.0, 0);

            var data = new double[hosking.OutputLength];
            hosking.SampleNoise(data);
        }
    }
}
