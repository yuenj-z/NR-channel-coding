import java.util.Random;

/**
 * The CommFunc class includes a few common functions 
 * that is needed in this project
 * 
 * @author Haiyue Zhang
 */

public class CommFunc {
    Random random = new Random();
    
    public double[] add(double[] first, double[] second){
        int length = first.length < second.length ? first.length : second.length;
        double[] result = new double[length];

        for (int i = 0; i < length; i++){
            result[i] = first[i] + second[i];
        }
        return result;
    }

    public double[] sub(double[] first, double[] second){
        int length = first.length < second.length ? first.length : second.length;
        double[] result = new double[length];

        for (int i = 0; i < length; i++){
            result[i] = first[i] - second[i];
        }
        return result;
    }

    public double[] create_msg(int k){
        double[] msg = new double[k];
        for (int i = 0; i < k; i++) {
            msg[i]  = random.nextBoolean() ? 1.0 : 0.0; 
        }
        return msg;
    }

    public double get_sign(double input){
        if (input > 0){
            return 1.0;
        }
        else {
            return -1.0;
        } 
    }

    public boolean contains(int[] array, int key) {
        for (int i : array) {
            if (i == key) {
                return true;
            }
        }
        return false;
    }
}
