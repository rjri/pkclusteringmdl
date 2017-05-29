package pkcluster;

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import java.io.FileWriter;
import java.io.IOException;
import javax.swing.UIManager;

public class Results {
 
	public JFrame frame;
	public JLabel lblTeseDoElson;
	public Output o;
	public double time;
	
	private JTextField textField_3;
	private JTextField textField_4;
	private JLabel lblQValue;
	private JTextField textField_5;
	private JTextField textField_6;
	private JTextField textField_7;
	private JTextField textField_8;
	private JTextArea textArea;


	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					Output o1=null;
					double t=0;
					Results window = new Results(o1,t);
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

	
	public Results(Output o1, double t) {
		o=o1;
		time=t;
		initialize();
	}
	
	public static void generateCsvFile(String sFileName,Output O)
	   {
		try
		{
		    FileWriter writer = new FileWriter(sFileName);
		    int i,l;
		    String Sw="",SS="",SC="",Si="";
		    writer.append("Results\n");
		    writer.append("Q value: "+O.Q+"\n");
		    writer.append("Clusters: "+O.cl+"\n");
		    
		    for(l=0;l<O.cl-1;l++){
		    	Sw=Sw+O.peso[l]+"; ";
		    	SS=SS+O.Sigma[l]+"; ";
		    	SC=SC+O.C[l][0]+"; "+O.C[l][1]+"; "+O.C[l][2]+"\n";
		    }
		    l=O.cl-1;
		    Sw=Sw+O.peso[l];
	    	SS=SS+O.Sigma[l];
	    	SC=SC+O.C[l][0]+"; "+O.C[l][1]+"; "+O.C[l][2]+"\n";
	    	
	    	writer.append("Weight: "+Sw+"\n");
	    	writer.append("Variance: "+SS+"\n");
	    	writer.append("\nParameters of clusters\n a; b1; b2:\n"+SC+"\n");
		    
	    	for(i=0;i<O.cluster.length-1;i++){
	    		Si=Si+O.cluster[i]+"; ";
	    	}
	    	i=O.cluster.length-1;
	    	Si=Si+O.cluster[i];
	    	writer.append("Clusters of individuals: "+Si+"\n");
	    	
	    	
	    	writer.flush();
		    writer.close();
		}
		catch(IOException e)
		{
		     System.out.println("ERROR");
		} 
	    }
	
	private void initialize() {
		frame = new JFrame();
		frame.setBounds(100, 100, 450, 350);
		frame.setForeground(new Color(255, 255, 255));
		frame.setFont(new Font("Dialog", Font.BOLD, 15));
		frame.setTitle("Unsupervised Learning Algorithm for Pharmacokinetic Responses\n");
		frame.setBackground(new Color(0, 0, 0));
		frame.getContentPane().setLayout(null);
		
		lblTeseDoElson = new JLabel("RESULTS");
		lblTeseDoElson.setBounds(9, 0, 210, 62);
		lblTeseDoElson.setForeground(new Color(0, 0, 153));
		lblTeseDoElson.setFont(new Font("Dialog", Font.BOLD, 26));
		frame.getContentPane().add(lblTeseDoElson);
		
		
		
		JButton btnSairEm = new JButton("Exit");
		btnSairEm.setBounds(330, 17, 98, 25);
		btnSairEm.setBackground(UIManager.getColor("OptionPane.errorDialog.titlePane.shadow"));
		btnSairEm.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				frame.setVisible(false);
			}
		});
		frame.getContentPane().add(btnSairEm);
		
			

		JButton btnEmRun = new JButton("Save");
		btnEmRun.setBounds(226, 17, 95, 25);
		btnEmRun.setBackground(UIManager.getColor("OptionPane.errorDialog.titlePane.shadow"));
		btnEmRun.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				generateCsvFile("Result.txt",o);
			}
		});
		frame.getContentPane().add(btnEmRun);
		
//		textField_1 = new JTextField();
//		textField_1.setBounds(425, 65, 93, 19);
//		frame.getContentPane().add(textField_1);
//		textField_1.setColumns(10);
		
		JLabel lblTempo = new JLabel("Time");
		lblTempo.setBounds(9, 47, 47, 15);
		frame.getContentPane().add(lblTempo);
		
		textField_3 = new JTextField();
		textField_3.setBounds(65, 45, 128, 19);
		frame.getContentPane().add(textField_3);
		textField_3.setColumns(10);
		
		textField_4 = new JTextField();
		textField_4.setBounds(9, 292, 419, 19);
		frame.getContentPane().add(textField_4);
		textField_4.setColumns(10);
		
		lblQValue = new JLabel("MDL");
		lblQValue.setBounds(202, 65, 53, 24);
		frame.getContentPane().add(lblQValue);
		
		textField_5 = new JTextField();
		textField_5.setBounds(262, 65, 166, 19);
		frame.getContentPane().add(textField_5);
		textField_5.setColumns(10);
		
		JLabel lblNClusters = new JLabel("Nb clusters");
		lblNClusters.setBounds(9, 65, 78, 24);
		frame.getContentPane().add(lblNClusters);
		
		textField_6 = new JTextField();
		textField_6.setBounds(90, 69, 103, 19);
		frame.getContentPane().add(textField_6);
		textField_6.setColumns(10);
		
		JLabel lblPesos = new JLabel("Weights");
		lblPesos.setBounds(9, 96, 74, 19);
		frame.getContentPane().add(lblPesos);
		
		textField_7 = new JTextField();
		textField_7.setBounds(90, 96, 338, 19);
		frame.getContentPane().add(textField_7);
		textField_7.setColumns(10);
		
		textField_8 = new JTextField();
		textField_8.setBounds(90, 122, 338, 19);
		frame.getContentPane().add(textField_8);
		textField_8.setColumns(10);
		
		JLabel lblSigma = new JLabel("Variance");
		lblSigma.setBounds(9, 122, 74, 19);
		frame.getContentPane().add(lblSigma);
		
		JLabel lblDadosDosPacientes = new JLabel("Subjects cluster");
		lblDadosDosPacientes.setBounds(9, 270, 210, 19);
		frame.getContentPane().add(lblDadosDosPacientes);
		
		
		textArea = new JTextArea();
		textArea.setBounds(9, 152, 419, 115);
		frame.getContentPane().add(textArea);
		
		textField_3.setText(time+" s");
		
		int i,l;
	    String Sw="",SS="",SC="",Si="";
	      
	    for(l=0;l<o.cl-1;l++){
	    	Sw=Sw+o.peso[l]+"; ";
	    	SS=SS+o.Sigma[l]+"; ";
	    	SC=SC+o.C[l][0]+"; "+o.C[l][1]+"; "+o.C[l][2]+"\n";
	    }
	    l=o.cl-1;
	    Sw=Sw+o.peso[l];
    	SS=SS+o.Sigma[l];
    	SC=SC+o.C[l][0]+"; "+o.C[l][1]+"; "+o.C[l][2]+"\n";
    	
    	for(i=0;i<o.cluster.length-1;i++){
    		Si=Si+o.cluster[i]+"; ";
    	}
    	i=o.cluster.length-1;
    	Si=Si+o.cluster[i];
		
		
		textArea.setText(SC);
		textField_7.setText(Sw);
		textField_8.setText(SS);
		textField_6.setText(o.cl+"");
		textField_5.setText(o.Q+"");
		textField_4.setText(Si);
		
		
	}

}
