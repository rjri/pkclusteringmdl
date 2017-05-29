package pkcluster;

import java.io.*;
public class Data{

	public Measurement top;
	public Measurement first;
	public int n;
	public double[] T;
	public int N;
	public double[][] M;
	public int[]ID;
	public double max;
	
	public Data(){ //Construtor da amostra.
		top=null;
		first=null;
		n=0;
		N=0;
		T=null;
		M=null;
		ID=null;
		max=0;
	}
	
	public void read(String args) throws Exception {
		//LE UM FICHEIRO E CRIA A AMOSTRA, FAZENDO ADD DE CADA LINHA DO FICHEIRO
		BufferedReader CSVfile = new BufferedReader(new FileReader(args));
		String dataRow=CSVfile.readLine();
		dataRow=CSVfile.readLine();
		while(dataRow!=null){
			String [] dataArray=dataRow.split(";");
			double[] intArray= new double[dataArray.length];
			int i=0;
			while(i<dataArray.length){
				intArray[i]=Double.parseDouble(dataArray[i]);
				i++;
			}
			add(intArray);
			dataRow=CSVfile.readLine();
		}
		CSVfile.close();
		end();
	}
	
	public double area(){
		double A=0;
		double m=max;
		A=m*T[T.length-1];
		return A;
	}
	
	public void actT(double t){
		int i,j,bool;
		double[] Ta;
		
		if(T==null){
			n++;
			T=new double[n];
			T[0]=t;}
		
		else{
			bool=0;
			for(i=0;i<T.length&&bool==0;i++){if(T[i]==t)bool=1;}
			
			if(bool==0){
				Ta=T;
				n++;
				T=new double[n];
				i=0;j=0;
				while(j<n){
					if(bool==0){
						if(i<n-1&&Ta[i]<t){T[j]=Ta[i];i++;}
						else{bool=1;T[j]=t;}
					}else{T[j]=Ta[i];i++;}
					j++;
				}
			}
			
		}
		
	}
	
	public void actID(int id){
		int i,j,bool;
		int[] I;
		
		if(ID==null){
			N++;
			ID=new int[N];
			ID[0]=id;}
		
		else{
			bool=0;
			for(i=0;i<ID.length&&bool==0;i++){if(ID[i]==id)bool=1;}
			
			if(bool==0){
				I=ID;
				N++;
				ID=new int[N];
				i=0;j=0;
				while(j<N){
					if(bool==0){
						if(i<N-1&&I[i]<id){ID[j]=I[i];i++;}
						else{bool=1;ID[j]=id;}
					}else{ID[j]=I[i];i++;}
					j++;
				}
			}
			
		}
		
	}
	
	public void add(int id, double t0, double y){ 
		Measurement aux = new Measurement(id,t0,y);
		int i;
		if(top!=null){
			top.next=aux;
			top=aux;}
		else{top=aux;
			first=aux;}

		actT(aux.data[1]);
		i=(int)aux.data[0];
		actID(i);
		if(y>max)max=y;
	}
	
	public void add(double[] x){ 
		Measurement aux = new Measurement(x);
		int i;
		double y=x[2];
		if(top!=null){
			top.next=aux;
			top=aux;}
		else{top=aux;
			first=aux;}

		actT(aux.data[1]);
		i=(int)aux.data[0];
		actID(i);
		if(y>max)max=y;
	}
		
	public int length(){
		return N;
	}
		
	public void join(Data B){
		if(top==null){
			first=B.first;
			top=B.top;}
		else{top.next=B.first;
			top=B.top;}
		}	
	
	public double[] element(int x){
		int i;
		Measurement aux;
		aux=first;
		for(i=0;i<x;i++)aux=aux.next;
		return aux.data;
	}
	
	public Measurement indice(int id){
		Measurement LA=null;
		Measurement aux=first;
		Measurement bux;
		while(aux!=null){
			if(aux.data[0]==id){
				bux=new Measurement(aux.data);
				bux.next=LA;
				LA=bux;}
			aux=aux.next;}
		return LA;
	}
	
	public void end(){
		M=new double[N][n];
		int id,i;
		double t,y;
		Measurement aux=first;
		while(aux!=null){
			y=aux.data[2];
			t=aux.data[1];
			id=(int)aux.data[0];
			i=0;while(T[i]!=t){i++;}
			M[id][i]=y;
			aux=aux.next;
		}
	}
	
	
}