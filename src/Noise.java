import java.util.Random;
import java.lang.Math;

/**
 * The Noise class generates random noise 
 * in the form of doubles 1.0 and 0.0
 * 
 * @author Haiyue Zhang
 */

public class Noise {
    public double[] generate_rand_num(int n, double mean, double variance){
        Random random = new Random();
        double[] noise = new double[n];

        for (int i = 0; i < n; i++){
            double a = random.nextDouble();
            double b = random.nextDouble();
            
            double rand = Math.sqrt(-2 * Math.log(a)) * Math.cos(2 * Math.PI * b);
            noise[i] = mean + Math.sqrt(variance) * rand;
        }
        return noise;
    }

    public int limit(double[] noise, double lower, double upper){
        int counter = 0;
        int n = noise.length;
        for (int i = 0; i < n; i++){
            if (noise[i] <= upper && noise[i] >= lower){
                counter++;
            }
        }
        return counter;
    }
}



