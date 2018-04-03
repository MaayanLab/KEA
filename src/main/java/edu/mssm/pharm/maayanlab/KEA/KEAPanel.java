package edu.mssm.pharm.maayanlab.KEA;

import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import edu.mssm.pharm.maayanlab.common.core.FileUtils;
import edu.mssm.pharm.maayanlab.common.core.SettingsChanger;
import edu.mssm.pharm.maayanlab.common.swing.FileDrop;
import edu.mssm.pharm.maayanlab.common.swing.UIUtils;

public class KEAPanel extends JPanel {

	private static final long serialVersionUID = -9052895312752852837L;
	
	static Logger log = Logger.getLogger(KEAPanel.class.getSimpleName());
	
	// JPanel
	private JPanel panel;

	// UI elements
	private JFileChooser openChooser, saveChooser;
	private JTextField openPath, savePath;
	private JTextArea inputTextArea, outputTextArea;
	private JButton openButton, runButton;
	private JComboBox interactionsCombo, sortByCombo, miningLevelCombo;
	private JTextField selectTopText;

	// Output
	private String output;
	
	public static void main(String[] args) {

		if (args.length == 0) {			
			// Schedule a job for the EDT
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					createAndShowGUI();
				}
			});
		}
		else {
			KEA.main(args);
		}
	}
	
	private static void createAndShowGUI() {
		// Try to use Nimbus look and feel
		try {            
            UIManager.setLookAndFeel("com.sun.java.swing.plaf.nimbus.NimbusLookAndFeel");
        } catch (Exception e) {
           log.warning("Nimbus: " + e);
        }
        
        // Create and set up the window
        JFrame appFrame = new JFrame("KEA - Kinase Enrichment Analysis");
        appFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        // Add content to the window
        KEAPanel appPanel = new KEAPanel();
        appFrame.setContentPane(appPanel);
        
        // Display the window
        appFrame.setResizable(false);
        appFrame.pack();
        appFrame.setVisible(true);
	}
	
	public KEAPanel() {
		this.setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		
		if (!Boolean.getBoolean("verbose"))
            log.setLevel(Level.WARNING);
		
		// Configure panel
		panel = this;

		// File choosers
		openChooser = new JFileChooser(System.getProperty("user.dir"));
		openChooser.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				File file = openChooser.getSelectedFile();
				if (file.canRead() && e.getActionCommand().equals(JFileChooser.APPROVE_SELECTION))
					setupIO(file);
			}
		});
		saveChooser = new JFileChooser(System.getProperty("user.dir"));
		saveChooser.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				File file = saveChooser.getSelectedFile();
				if (file != null && e.getActionCommand().equals(JFileChooser.APPROVE_SELECTION)) {
					if (!file.getName().endsWith(".csv")) {
						file = new File(file.getAbsolutePath() + ".csv");
						saveChooser.setSelectedFile(file);
					}
					
					savePath.setText(file.getAbsolutePath());
				}
			}
		});
		
		// Select input file button
		JButton openFileButton = new JButton("Input Proteins");
		openFileButton.setPreferredSize(new Dimension(300, 30));
		openFileButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				openChooser.showOpenDialog(panel);
			}
		});
		JButton saveFileButton = new JButton("Output Kinases");
		saveFileButton.setPreferredSize(new Dimension(300, 30));
		saveFileButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				saveChooser.showSaveDialog(panel);
			}
		});
		
		// Text Fields
		openPath = new JTextField();
		savePath = new JTextField();
		
		// File Drop
		new FileDrop(openPath, new FileDrop.Listener() {
			public void filesDropped(File[] files) {
				if (files[0].canRead()) {
					setupIO(files[0]);
					openChooser.setSelectedFile(files[0]);
				}
			}
		});
		
		// Scroll panes
		inputTextArea = new JTextArea(20, 20);
		JScrollPane inputTextPane = new JScrollPane(inputTextArea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		outputTextArea = new JTextArea(20, 20);
		JScrollPane outputTextPane = new JScrollPane(outputTextArea, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
		
		// File Drop
		new FileDrop(inputTextArea, new FileDrop.Listener() {
			public void filesDropped(File[] files) {
				if (files[0].canRead()) {
					setupIO(files[0]);
					openChooser.setSelectedFile(files[0]);
				}
			}
		});
		
		// Open results
		openButton = new JButton("View Results");
		openButton.setEnabled(false);
		openButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				try {
					Desktop.getDesktop().open(new File(output));
				} catch (Exception e1) {
					JOptionPane.showMessageDialog(panel, "Unable to open " + output, "Unable to open file", JOptionPane.ERROR_MESSAGE);
				}
			}
		});
		
		// Start button
		runButton = new JButton("Find Kinases");
		runButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				output = savePath.getText();
				ArrayList<String> inputList = UIUtils.getTextAreaText(inputTextArea);
				
				try {
					if (!output.equals("") && FileUtils.validateList(inputList)) {
						KEA kea = new KEA();
						
						setSettings(kea);
						kea.run(inputList);
						UIUtils.setTextAreaText(outputTextArea, kea.getTopRankedList(Integer.parseInt(selectTopText.getText())));
						kea.writeFile(output);
						enableOutput(output);
					}
					else {
						JOptionPane.showMessageDialog(panel, "No save location specified.", "No Save Location", JOptionPane.WARNING_MESSAGE);
					}
				} catch (ParseException e1) {
					if (e1.getErrorOffset() == -1)
						JOptionPane.showMessageDialog(panel, "Input list is empty.", "Invalid Input", JOptionPane.WARNING_MESSAGE);
					else
						JOptionPane.showMessageDialog(panel, e1.getMessage() + " at line " + (e1.getErrorOffset() + 1) +" is not a valid Entrez Gene Symbol.", "Invalid Input", JOptionPane.WARNING_MESSAGE);
				}
			}
		});
		
		// Advanced Settings
		JLabel interactionsLabel1 = new JLabel("Include");
		String[] interactionsOptions = {"kinase-protein interactions", "phosphorylation reactions", "kinase-protein and phosphorylation interactions"};
		interactionsCombo = new JComboBox(interactionsOptions);
		interactionsCombo.setSelectedIndex(2);
		JLabel interactionsLabel2 = new JLabel("in the background database");
		JPanel interactionsBox = new JPanel();
		interactionsBox.add(interactionsLabel1);
		interactionsBox.add(interactionsCombo);
		interactionsBox.add(interactionsLabel2);
		
		JLabel sortByLabel = new JLabel("Sort by");
		String[] sortingOptions = {"p-value", "rank", "combined score"};
		sortByCombo = new JComboBox(sortingOptions);
		sortByCombo.setSelectedIndex(2);
		JPanel sortByBox = new JPanel();
		sortByBox.add(sortByLabel);
		sortByBox.add(sortByCombo);
		
		JLabel selectTopLabel1 = new JLabel("Select top");
		selectTopText = new JTextField("10");
		JLabel selectTopLabel2 = new JLabel("kinases");
		JPanel selectTopBox = new JPanel();
		selectTopBox.add(selectTopLabel1);
		selectTopBox.add(selectTopText);
		selectTopBox.add(selectTopLabel2);
		
		JLabel miningLevelLabel = new JLabel("Resolve down to");
		String[] choice = {"kinases", "kinase families", "kinase classes"};
		miningLevelCombo = new JComboBox(choice);
		miningLevelCombo.setSelectedIndex(0);
		JPanel miningLevelBox = new JPanel();
		miningLevelBox.add(miningLevelLabel);
		miningLevelBox.add(miningLevelCombo);
		
		// Input and output box
		JPanel ioBox = new JPanel();
		ioBox.setLayout(new GridLayout(2,2));
		ioBox.add(openFileButton);
		ioBox.add(saveFileButton);
		ioBox.add(openPath);	
		ioBox.add(savePath);
		
		// Panes
		JPanel textPanesBox = new JPanel();
		textPanesBox.setLayout(new BoxLayout(textPanesBox, BoxLayout.LINE_AXIS));
		textPanesBox.add(inputTextPane);
		textPanesBox.add(outputTextPane);
		
		// Button box
		JPanel buttonBox = new JPanel();
		buttonBox.setLayout(new BoxLayout(buttonBox, BoxLayout.LINE_AXIS));
		buttonBox.add(runButton);
		buttonBox.add(openButton);
		
		// Advanced settings box
		JPanel advancedSettingsBox = new JPanel();
		advancedSettingsBox.setLayout(new BoxLayout(advancedSettingsBox, BoxLayout.PAGE_AXIS));
		advancedSettingsBox.setBorder(BorderFactory.createTitledBorder("Advanced Settings"));
		advancedSettingsBox.add(interactionsBox);
		advancedSettingsBox.add(sortByBox);
		advancedSettingsBox.add(selectTopBox);
		advancedSettingsBox.add(miningLevelBox);
		
		this.add(ioBox);
		this.add(textPanesBox);
		this.add(Box.createRigidArea(new Dimension(0,10)));
		this.add(buttonBox);
		this.add(advancedSettingsBox);
	}
	
	public void setSettings(SettingsChanger changer) {
		switch (interactionsCombo.getSelectedIndex()) {
		case 0: changer.setSetting(KEA.KINASE_INTERACTIONS, KEA.KINASE_PROTEIN); break;
		case 1: changer.setSetting(KEA.KINASE_INTERACTIONS, KEA.PHOSPHORYLATION); break;
		case 2: changer.setSetting(KEA.KINASE_INTERACTIONS, KEA.BOTH_TYPES); break;
		}
		switch (sortByCombo.getSelectedIndex()) {
		case 0: changer.setSetting(KEA.SORT_BY, KEA.PVALUE); break;
		case 1: changer.setSetting(KEA.SORT_BY, KEA.RANK); break;
		case 2: changer.setSetting(KEA.SORT_BY, KEA.COMBINED_SCORE); break;
		}
		switch (miningLevelCombo.getSelectedIndex()) {
		case 0: changer.setSetting(KEA.RESOLUTION_LEVEL, KEA.KINASE_LEVEL); break;
		case 1: changer.setSetting(KEA.RESOLUTION_LEVEL, KEA.KINASE_FAMILY_LEVEL); break;
		case 2: changer.setSetting(KEA.RESOLUTION_LEVEL, KEA.KINASE_GROUP_LEVEL); break;
		}
		changer.setSetting("number_of_top_kinases", selectTopText.getText());
	}
	
	public void setInputTextArea(Collection<String> list) {
		UIUtils.setTextAreaText(inputTextArea, list);
	}
	
	public void setOutputTextArea(Collection<String> list) {
		UIUtils.setTextAreaText(outputTextArea, list);
	}
	
	private void setupIO(File inputFile) {
		openPath.setText(inputFile.getAbsolutePath());
		UIUtils.setTextAreaText(inputTextArea, FileUtils.readFile(inputFile));
		
		File outputFile = new File(System.getProperty("user.dir"), FileUtils.stripFileExtension(inputFile.getName()) + ".results_kinase.csv");
		saveChooser.setSelectedFile(outputFile);
		savePath.setText(outputFile.getAbsolutePath());
	}
	
	public void enableOutput(String output) {
		savePath.setText(output);
		this.output = output;
		if (Desktop.isDesktopSupported() && Desktop.getDesktop().isSupported(Desktop.Action.OPEN))
			openButton.setEnabled(true);
	}
}