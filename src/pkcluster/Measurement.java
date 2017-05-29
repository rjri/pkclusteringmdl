package pkcluster;

public class Measurement {
	public double[] data;
	public Measurement next;
	
	public Measurement (int id, double time, double concentration){
		next=null;
		data=new double[3];
		data[0]=id;
		data[1]=time;
		data[2]=concentration;
	}
	
	public Measurement (double[] x){
		next=null;
		data=x;
	}
}
