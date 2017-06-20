package pkcluster;

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Font;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JButton;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.SwingConstants;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;


public class EMRun {

	public JFrame frame;
	public JLabel lblTeseDoElson;
	private JSlider slider;
	private JTextField textField_1;
	public int Max;
	public int n;
	public Data A;
	public Output o;
	public double time;
	
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					int Max1=15;
					int n1=10;
					Data A1=new Data();
					EMRun window = new EMRun(A1, Max1, n1);
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	public EMRun(Data A1, int Max1, int n1) {
		A=A1;Max=Max1;n=n1;
		initialize();
	}

	private void initialize() {
		frame = new JFrame();
		frame.setBounds(100, 100, 620, 180);
		frame.setForeground(new Color(50, 50, 50));
		frame.setFont(new Font("Dialog", Font.BOLD, 15));
		frame.setTitle("Unsupervised Learning Algorithm for Pharmacokinetic Responses\n");
		frame.setBackground(new Color(255, 255, 255));
		frame.getContentPane().setLayout(null);
		
		lblTeseDoElson = new JLabel("Enter EM parameters");
		lblTeseDoElson.setBounds(9, 0, 509, 50);
		lblTeseDoElson.setForeground(new Color(0, 0, 153));
		lblTeseDoElson.setFont(new Font("Dialog", Font.BOLD, 26));
		frame.getContentPane().add(lblTeseDoElson);
		
		final JTextPane textPane = new JTextPane();
		textPane.setBounds(9, 92, 597, 47);
		textPane.setBackground(Color.GRAY);
		textPane.setForeground(Color.BLUE);
		textPane.setFont(new Font("Dialog", Font.BOLD, 24));
		frame.getContentPane().add(textPane);
		
		JButton btnSairEm = new JButton("Exit");
		btnSairEm.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				frame.setVisible(false);
			}
		});
		btnSairEm.setBounds(514, 17, 92, 25);
		frame.getContentPane().add(btnSairEm);
		
		slider = new JSlider();
		slider.setBounds(237, 44, 273, 50);
		slider.setValue(10);
		slider.setOrientation(SwingConstants.HORIZONTAL);
		slider.setPaintLabels(true);
		slider.setMaximum(Max);
		slider.setMinimum(n);
		frame.getContentPane().add(slider);
		slider.setFont(new Font("Tahoma",Font.BOLD,12));
        slider.setMajorTickSpacing((int)(Max-n)/5);
        slider.setMinorTickSpacing((int)(Max-n)/100);
        slider.setPaintLabels(true);
        slider.setPaintTicks(true);
        slider.setPaintTrack(true);
        slider.setAutoscrolls(true);
		
        JLabel lblNClusters = new JLabel("Num clusters");
		lblNClusters.setBounds(159, 44, 93, 22);
		frame.getContentPane().add(lblNClusters);
		
		JLabel lblNInputs = new JLabel("Rand Init");
		lblNInputs.setBounds(9, 48, 70, 15);
		frame.getContentPane().add(lblNInputs);
		
		textField_1 = new JTextField();
		textField_1.setBounds(81, 46, 57, 19);
		frame.getContentPane().add(textField_1);
		textField_1.setColumns(10);
		
		JButton btnEmRun = new JButton("Run");
		btnEmRun.setBounds(514, 44, 92, 25);
		btnEmRun.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				
				
				Output o1,o2;
				int k,K;
				
				double Q=-150;
				long startTime,endTime;
				double duration;
				final int M0;
				try{ 
					K=Integer.parseInt(textField_1.getText());
					startTime = System.currentTimeMillis();
					//M0= 3;/* ASMC slider.getValue();*/
					M0=slider.getValue();
					EMAlgorithm E;
					o2=null;
					/*for(k=0;k<K;k++){
						E=new EMAlgorithm(A.T,A.M,M0);
						E.inicialize(A);
						o1=E.runEM();
						if(k==0||o1.Q>Q){Q=o1.Q;o2=o1;}
						
						//HERE : asmc
						int i;
					    String Si="";
						for(i=0;i<o2.cluster.length-1;i++){
				    		Si=Si+o2.cluster[i]+"; ";
				    	}
				    	i=o2.cluster.length-1;
				    	Si=Si+o2.cluster[i];
						System.out.println("M: "+ EMAlgorithm.maxMDLnumclusttotal+">> "+Si);
						
						
					}*/
					
					ExecutorService executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
					List<Future<Output>> futures = new ArrayList<Future<Output>>();
					for(k=0;k<K;k++){
						Callable<Output> callable = new Callable<Output>(){
							public Output call(){
								EMAlgorithm E=new EMAlgorithm(A.T,A.M,M0);
								E.inicialize(A);
								Output op=E.runEM();
								return op;
							}
						};
						futures.add(executorService.submit(callable));
					}
					executorService.shutdown();
					executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
					int num=0;
					for(Future<Output> future : futures){
						o1=future.get();
						System.out.println("num: "+num);o1.imprime();
						if(num==0||o1.Q>Q){Q=o1.Q;o2=o1;System.out.println("^ new best");}
						num++;
					}
					
					//System.out.println("Max MDL result: "+EMAlgorithm.maxMDLtotal+"; Number of Clusters: "+EMAlgorithm.maxMDLnumclusttotal);
					//o=new Output(EMAlgorithm.maxMDL_C_params,EMAlgorithm.maxMDLtotal,k,EMAlgorithm.maxMDL_Sigma,EMAlgorithm.maxMDL_w,EMAlgorithm.maxMDLnumclusttotal,EMAlgorithm.maxMDL_cluster);
					o2.imprime();
					
					endTime = System.currentTimeMillis();
					duration = endTime - startTime;
					duration=duration*0.001/k;
					o=o2;
					time=duration;
					textPane.setText("DONE!");
					textPane.setForeground(Color.WHITE);
					
					Results window = new Results(o,time);
					window.frame.setVisible(true);
					
				} catch (Exception e1){
					textPane.setText("Errors!");
					textPane.setForeground(Color.RED);
					};
			}
		});
		frame.getContentPane().add(btnEmRun);
		
	}
}
