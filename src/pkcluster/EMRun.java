package pkcluster;

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Font;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JButton;
import javax.swing.JTextField;
import javax.swing.JTextPane;
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
	public JLabel lbl;
	private JTextField textField_1, textField_2, textField_3;
	public int Max;
	public int n;
	public Data A;
	public Output o;
	public double time;
	public boolean nml;
	
	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					int Max1=15;
					int n1=10;
					Data A1=new Data();
					boolean nml1=false;
					EMRun window = new EMRun(A1, Max1, n1,nml1);
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}
	
	public EMRun(Data A1, int Max1, int n1, boolean nml1) {
		A=A1;Max=Max1;n=n1;nml=nml1;
		initialize();
	}

	private void initialize() {
		frame = new JFrame();
		frame.setBounds(100, 100, 620, 225);
		frame.setForeground(new Color(50, 50, 50));
		frame.setFont(new Font("Dialog", Font.BOLD, 15));
		frame.setTitle("Unsupervised Learning Algorithm for Pharmacokinetic Responses\n");
		frame.setBackground(new Color(255, 255, 255));
		frame.getContentPane().setLayout(null);
		
		lbl = new JLabel("Enter clustering parameters");
		lbl.setBounds(9, 0, 509, 50);
		lbl.setForeground(new Color(0, 0, 153));
		lbl.setFont(new Font("Dialog", Font.BOLD, 26));
		frame.getContentPane().add(lbl);
		
		final JTextPane textPane = new JTextPane();
		textPane.setBounds(9, 137, 597, 47);
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
		
        JLabel lblNClusters = new JLabel("Maximum number of clusters:");
		lblNClusters.setBounds(9, 98, 210, 15);
		frame.getContentPane().add(lblNClusters);
		
		JLabel lblNClusters2 = new JLabel("Minimum number of clusters:");
		lblNClusters2.setBounds(9, 73, 210, 15);
		frame.getContentPane().add(lblNClusters2);
		
		JLabel lblNInputs = new JLabel("Random initializations:");
		lblNInputs.setBounds(9, 48, 170, 15);
		frame.getContentPane().add(lblNInputs);
		
		textField_1 = new JTextField();
		textField_1.setBounds(230, 46, 57, 19);
		frame.getContentPane().add(textField_1);
		textField_1.setColumns(10);
		
		textField_2 = new JTextField();
		textField_2.setBounds(230, 96, 57, 19);
		frame.getContentPane().add(textField_2);
		textField_2.setColumns(10);
		
		textField_3 = new JTextField();
		textField_3.setBounds(230, 71, 57, 19);
		frame.getContentPane().add(textField_3);
		textField_3.setColumns(10);
		
		JButton btnEmRun = new JButton("Run");
		btnEmRun.setBounds(514, 44, 92, 25);
		btnEmRun.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				
				
				Output o1,o2;
				int k,K;
				
				double Q=-150;
				long startTime,endTime;
				double duration;
				final int M0, m0;
				try{ 
					K=Integer.parseInt(textField_1.getText());
					startTime = System.currentTimeMillis();
					M0=Integer.parseInt(textField_2.getText());
					m0=Integer.parseInt(textField_3.getText());
					if(m0>M0) throw new Exception();
					o2=null;
					
					ExecutorService executorService = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
					List<Future<Output>> futures = new ArrayList<Future<Output>>();
					for(k=0;k<K;k++){
						Callable<Output> callable = new Callable<Output>(){
							public Output call(){
								EMAlgorithm E=new EMAlgorithm(A.T,A.M,M0,nml,m0);
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
						if(num==0||o1.Q>Q){Q=o1.Q;o2=o1;}
						num++;
					}
					
					o2.imprime();
					
					endTime = System.currentTimeMillis();
					duration = endTime - startTime;
					duration=duration*0.001;
					o=o2;
					time=duration;
					textPane.setText("DONE!");
					textPane.setForeground(Color.WHITE);
					
					Results window = new Results(o,time,nml);
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
