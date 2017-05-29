package pkcluster;

public class Output {
	public double[][] C;
	public double Q;
	public int passos;
	public double[] Sigma;
	public double[] peso;
	public int cl; 
	public Output next;
	int[] cluster;
	
	public Output(double[][] C0, double Q0, int p, double[] S,double[] w, int c,int[] cluster0){
		C=C0;
		Q=Q0;
		passos=p;
		Sigma=S;
		peso=w;
		cl=c;
		next=null;
		cluster=cluster0;
		
	}
	
	
	public void imprime(){
		int l;
		System.out.println("OUTPUT: "+passos);
		System.out.println();
		
		System.out.println("Sigma: ");
		for(l=0;l<cl;l++)System.out.print(Sigma[l]+"; ");
		System.out.println();
		
		System.out.println("Pesos: ");
		for(l=0;l<cl;l++)System.out.print(peso[l]+"; ");
		System.out.println();
		
		System.out.println("Number of clusters: "+cl);
		
		System.out.println("a; ke; ka");
		for(l=0;l<cl;l++)
			System.out.println("{"+C[l][0]+", "+C[l][1]+", "+C[l][2]+"},");		
		System.out.println("MDL value: "+Q);
	}
	
	public String string(){
		String S,S1;
		int l;
		
		S="";		
		S1="Sigma: {";
		for(l=0;l<cl;l++)S1=S1+Sigma[l]+", ";
		S1=S1+"}\n";
		S=S+S1;
		
		S1="Pesos: {";
		for(l=0;l<cl;l++)S1=S1+peso[l]+", ";
		S1=S1+"}\n";
		S=S+S1;
		
		S1="Number of clusters: "+cl+"\n";
		S=S+S1;
		
		S1="{";
		for(l=0;l<cl;l++)
			S1=S1+"{"+C[l][0]+", "+C[l][1]+", "+C[l][2]+"},";		
		S1=S1+"}\n";
		S=S+S1;
		S=S+"Q-value: "+Q+"\n";
		
		return S;
	}
	
}
