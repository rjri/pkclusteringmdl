package pkcluster;

import java.awt.Color;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.JTextPane;

public class Main {

	private JFrame frame;
	private JTextField textField;
	private JLabel lbl;
	int Max=15;
	int n=10;
	Data A=new Data();

	public static void main(String[] args) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
					Main window = new Main();
					window.frame.setVisible(true);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}


	public Main() {
		initialize();
	}


	private void initialize() {
		frame = new JFrame();
		frame.setBounds(100, 100, 600, 180);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().setLayout(null);
		frame.setForeground(new Color(255, 255, 255));
		frame.setFont(new Font("Dialog", Font.BOLD, 15));
		frame.setTitle("Unsupervised Learning Algorithm for Pharmacokinetic Responses\n");
		frame.setBackground(new Color(0, 0, 0));
		
		lbl = new JLabel("Insert PK data file");
		lbl.setBounds(9, 0, 509, 50);
		lbl.setForeground(new Color(0, 0, 153));
		lbl.setFont(new Font("Dialog", Font.BOLD, 26));
		frame.getContentPane().add(lbl);
		
		JLabel lblDadosDeAmostra = new JLabel("Data file in CSV format:");
		lblDadosDeAmostra.setBounds(9, 61, 170, 27);
		lblDadosDeAmostra.setBackground(new Color(255, 255, 255));
		frame.getContentPane().add(lblDadosDeAmostra);
		
		textField = new JTextField();
		textField.setBounds(183, 65, 153, 19);
		frame.getContentPane().add(textField);
		textField.setColumns(10);
		
		final JTextPane textPane = new JTextPane();
		textPane.setBounds(9, 92, 577, 47);
		textPane.setBackground(Color.BLACK);
		textPane.setForeground(Color.BLUE);
		textPane.setFont(new Font("Dialog", Font.BOLD, 24));
		frame.getContentPane().add(textPane);
		
		final JRadioButton mdlbutton = new JRadioButton("MDL");
		mdlbutton.setBounds(389, 30, 97, 20);
		mdlbutton.setSelected(true);
		frame.getContentPane().add(mdlbutton);
		
		final JRadioButton nmlbutton = new JRadioButton("NML");
		nmlbutton.setBounds(489, 30, 97, 20);
		frame.getContentPane().add(nmlbutton);
		
		ButtonGroup group = new ButtonGroup();
		group.add(mdlbutton);
		group.add(nmlbutton);
		
		JButton btnAmostra = new JButton("Read data");
		btnAmostra.setBounds(344, 65, 145, 19);
		btnAmostra.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				try{
					boolean nml;
					if(nmlbutton.isSelected()){
						nml=true;
					}else{
						nml=false;
					}
					
					A.read(textField.getText());
					Max=A.N;
					n=Max;n=n/3;
					if(Max<15)Max=15;
					if(n<10)n=10;
					
					textPane.setText("Data Read!");
					textPane.setForeground(Color.WHITE);
					EMRun window2 = new EMRun(A, Max, n, nml);
					window2.frame.setVisible(true);
					
					} catch (Exception e2){
					textPane.setText("File not found!");
					textPane.setForeground(Color.RED);
					};
			}
		});
		frame.getContentPane().add(btnAmostra);
		
		
		
		JButton btnClear = new JButton("Clear");
		btnClear.setBounds(501, 65, 85, 19);
		btnClear.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent arg0) {
				textField.setText("");
				}
		});
		frame.getContentPane().add(btnClear);
	}
}
