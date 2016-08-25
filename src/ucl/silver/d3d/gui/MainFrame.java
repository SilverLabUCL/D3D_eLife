package ucl.silver.d3d.gui;

import ucl.silver.d3d.core.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/**
 * <p>Title: D3D</p>
 * <p>Description: 3D Diffusion-Reaction Simulator</p>
 * <p>Copyright: Copyright (c) 2016</p>
 * <p>Company: The Silver Lab at University College London</p>
 * @author Jason S. Rothman
 * @version 1.6
 */

public class MainFrame
        extends JFrame {

    JPanel contentPane;

    JFileChooser fc = new JFileChooser();

    JMenuBar jManuMain = new JMenuBar();

    JMenu jMenuD3D = new JMenu();
    JMenuItem jMenuD3DAbout = new JMenuItem();
    JMenuItem jMenuD3DExit = new JMenuItem();

    JMenu jMenuProject = new JMenu();
    JMenuItem jMenuProjectInit = new JMenuItem();
    JMenuItem jMenuProject_FD_FRAP = new JMenuItem();
    JMenuItem jMenuProject_MC_FRAP = new JMenuItem();
    JMenuItem jMenuProject_MC_FRAP_Movie = new JMenuItem();
    JMenuItem jMenuProject_MC_MSD = new JMenuItem();
    JMenuItem jMenuProject_MC_AZ = new JMenuItem();

    JTabbedPane tabbedPane = new JTabbedPane();
    JPanel tabParam = new JPanel();
    JPanel tab2Dview = new JPanel();
    JPanel tabLog = new JPanel();

    public PanelParams panelParams;
    public Panel2D panel2D = null;
    public PanelLog panelLog;
    public JPanel panelCompute = new JPanel();

    JButton jButtonComputeJava = new JButton();
    JButton jButtonPause = new JButton();
    JButton jButtonPreview = new JButton();
    JButton jButtonCancel = new JButton();
    JButton jButtonInit = new JButton();

    FlowLayout flowLayoutCompute = new FlowLayout();
    BorderLayout borderLayoutParam = new BorderLayout();
    BorderLayout borderLayout2Dview = new BorderLayout();
    BorderLayout borderLayoutLog = new BorderLayout();
    BorderLayout borderLayoutContent = new BorderLayout();
    FlowLayout flowLayout1 = new FlowLayout();
    
    ActionListener projectActionListener;

    // construct the frame
    public MainFrame() {

        enableEvents(AWTEvent.WINDOW_EVENT_MASK);

        try {
            jbInit();
            validate();
        } catch (Exception e) {
            //e.printStackTrace();
        }

    }

    // Component initialization
    public final void jbInit() throws Exception {

        contentPane = (JPanel) this.getContentPane();
        contentPane.setLayout(borderLayoutContent);
        contentPane.setBorder(BorderFactory.createLoweredBevelBorder());
        setSize(new Dimension(1000, 800));
        //this.setTitle("D3D");

        jMenuD3D.setText("D3D");

        jMenuD3DAbout.setText("About");
        jMenuD3DAbout.addActionListener(new MainFrame_jMenuD3DAbout_ActionAdapter(this));
        jMenuD3DExit.setText("Quit");
        jMenuD3DExit.addActionListener(new MainFrame_jMenuD3DExit_ActionAdapter(this));
        
        projectActionListener = new MainFrame_jMenuProject_actionAdapter(this);

        jMenuProject.setText("eLife");
        jMenuProjectInit.setText("Init");
        jMenuProjectInit.addActionListener(projectActionListener);
        jMenuProject_FD_FRAP.setText("Init Finite Difference FRAP");
        jMenuProject_FD_FRAP.addActionListener(projectActionListener);
        jMenuProject_MC_FRAP.setText("Init Monte Carlo FRAP");
        jMenuProject_MC_FRAP.addActionListener(projectActionListener);
        jMenuProject_MC_FRAP_Movie.setText("Init Monte Carlo FRAP Movie");
        jMenuProject_MC_FRAP_Movie.addActionListener(projectActionListener);
        jMenuProject_MC_MSD.setText("Init Monte Carlo MSD");
        jMenuProject_MC_MSD.addActionListener(projectActionListener);
        jMenuProject_MC_AZ.setText("Init Monte Carlo AZ");
        jMenuProject_MC_AZ.addActionListener(projectActionListener);

        jButtonComputeJava.setPreferredSize(new Dimension(100, 25));
        jButtonComputeJava.setActionCommand("Compute");
        jButtonComputeJava.setText("Run");
        jButtonComputeJava.addActionListener(new MainFrame_jButtonComputeJava_actionAdapter(this));

        panelParams = new PanelParams();
        panel2D = new Panel2D();
        panelLog = new PanelLog();

        tabParam.setLayout(borderLayoutParam);

        jButtonPause.setPreferredSize(new Dimension(100, 25));
        jButtonPause.setActionCommand("Pause");
        jButtonPause.setText("Pause");
        jButtonPause.addActionListener(new MainFrame_jButtonPause_actionAdapter(this));

        jButtonPreview.setPreferredSize(new Dimension(100, 25));
        jButtonPreview.setActionCommand("Preview");
        jButtonPreview.setText("Preview");
        jButtonPreview.addActionListener(new MainFrame_jButtonPreview_actionAdapter(this));

        jButtonCancel.setAlignmentY((float) 0.5);
        jButtonCancel.setPreferredSize(new Dimension(100, 25));
        jButtonCancel.setActionCommand("Cancel");
        jButtonCancel.setText("Cancel");
        jButtonCancel.addActionListener(new MainFrame_jButtonCancel_actionAdapter(this));

        jButtonInit.setAlignmentY((float) 0.5);
        jButtonInit.setPreferredSize(new Dimension(100, 25));
        jButtonInit.setActionCommand("Init");
        jButtonInit.setText("Init");
        jButtonInit.addActionListener(new MainFrame_jButtonInit_actionAdapter(this));

        tabParam.add(panelParams, BorderLayout.CENTER);

        tab2Dview.setLayout(borderLayout2Dview);
        tab2Dview.add(panel2D, BorderLayout.CENTER);

        tabLog.setLayout(borderLayoutLog);
        tabLog.add(panelLog, BorderLayout.CENTER);

        tabbedPane.add(tabParam, "Parameters");
        tabbedPane.add(tab2Dview, "2D view");
        tabbedPane.add(tabLog, "Log");
        tabbedPane.addMouseListener(new MainFrame_tabbedPane_mouseAdapter(this));

        panelCompute.setLayout(flowLayoutCompute);
        panelCompute.add(jButtonInit, null);
        panelCompute.add(jButtonComputeJava, null);
        panelCompute.add(jButtonPause, null);
        panelCompute.add(jButtonPreview, null);
        panelCompute.add(jButtonCancel, null);
        contentPane.add(panelCompute, BorderLayout.SOUTH);

        jMenuD3D.add(jMenuD3DAbout);
        jMenuD3D.add(jMenuD3DExit);

        jManuMain.add(jMenuD3D);
        jManuMain.add(jMenuProject);

        setJMenuBar(jManuMain);
        contentPane.add(tabbedPane, BorderLayout.CENTER);

        initTabs();

    }

    public void updateMenuProject() {

        String projectList[] = Master.projectList();

        int numProjects = projectList.length;

        jMenuProject.removeAll();
        jMenuProject.add(jMenuProject_FD_FRAP);
        jMenuProject.add(jMenuProject_MC_FRAP);
        jMenuProject.add(jMenuProject_MC_FRAP_Movie);
        jMenuProject.add(jMenuProject_MC_MSD);
        jMenuProject.add(jMenuProject_MC_AZ);
        jMenuProject.add("---");
        jMenuProject.add(jMenuProjectInit);

    }

    public void initTabs() {

        if (panelParams != null) {
            panelParams.updateControls();
        }

        updateMenuProject();

        repaint();

    }

    public void tabbedPaneClicked(MouseEvent e) {
        
        switch(tabbedPane.getSelectedIndex()) {
            case 0:
                break;
            case 1:
                panel2D.updateControls();
                break;
            case 2:
                break;
        }

    }

    // File | Exit action performed
    public void jMenuD3DExit_actionPerformed(ActionEvent e) {
        System.exit(0);
    }

    // Help | About action performed
    public void jMenuD3DAbout_actionPerformed(ActionEvent e) {
        AboutD3D jDialog = new AboutD3D(this);
        Dimension dlgSize = jDialog.getPreferredSize();
        Dimension frmSize = getSize();
        Point loc = getLocation();
        jDialog.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x,
                (frmSize.height - dlgSize.height) / 2 + loc.y);
        jDialog.setModal(true);
        jDialog.pack();
        jDialog.setVisible(true);
    }

    // Overridden so we can exit when window is closed
    @Override
    protected void processWindowEvent(WindowEvent e) {
        super.processWindowEvent(e);
        if (e.getID() == WindowEvent.WINDOW_CLOSING) {
            jMenuD3DExit_actionPerformed(null);
        }
    }

    void jButtonComputeJava_actionPerformed(ActionEvent e) {
        Master.project.simulationStart(false);
    }

    void jButtonPause_actionPerformed(ActionEvent e) {
        Master.project.simulationPause();
    }

    void jButtonPreview_actionPerformed(ActionEvent e) {
        Master.project.simulationStart(true);
    }

    void jButtonCancel_actionPerformed(ActionEvent e) {
        Master.project.simulationCancel();
    }

    void jButtonInit_actionPerformed(ActionEvent e) {

        if (Master.project.monteCarlo != null) {
            Master.project.monteCarlo.initAll();
        }

        Master.project.init();

    }

    void jMenuProject_actionPerformed(ActionEvent e) {
        
        boolean error;
        String initClass = "", initFxn = "";

        String command = e.getActionCommand();
        
        if (command.equalsIgnoreCase("Init")) {
            error = Master.initProjectClassAndFunction("", "");
            return;
        }

        if (command.equalsIgnoreCase("Init Finite Difference FRAP")) {
            initClass = "InitFrapJason";
            initFxn = "init_FD_FRAP";
        }
        
        if (command.equalsIgnoreCase("Init Monte Carlo FRAP")) {
            initClass = "InitFrapJason";
            initFxn = "init_MC_FRAP_Prompt";
        }
        
        if (command.equalsIgnoreCase("Init Monte Carlo FRAP Movie")) {
            initClass = "InitFrapJason";
            initFxn = "init_MC_FRAP_Movie";
        }
        
        if (command.equalsIgnoreCase("Init Monte Carlo MSD")) {
            initClass = "InitFrapJason";
            initFxn = "init_MC_MSD_Prompt";
        }
        
        if (command.equalsIgnoreCase("Init Monte Carlo AZ")) {
            initClass = "InitAZEMJason";
            initFxn = "init_MC_AZ_Prompt";
        }
        
        Master.initProjectClass(initClass);
        Master.project.nullArrays();
        
        error = Master.initProjectClassAndFunction(initClass, initFxn);

    }
    
}

class MainFrame_jMenuD3DExit_ActionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuD3DExit_ActionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuD3DExit_actionPerformed(e);
    }
}

class MainFrame_jMenuD3DAbout_ActionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuD3DAbout_ActionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuD3DAbout_actionPerformed(e);
    }
}

class MainFrame_jButtonComputeJava_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonComputeJava_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonComputeJava_actionPerformed(e);
    }
}

class MainFrame_jMenuProject_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jMenuProject_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jMenuProject_actionPerformed(e);
    }
}

class MainFrame_tabbedPane_mouseAdapter
        extends java.awt.event.MouseAdapter {

    MainFrame adaptee;

    MainFrame_tabbedPane_mouseAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    @Override
    public void mouseClicked(MouseEvent e) {
        adaptee.tabbedPaneClicked(e);
    }

}

class MainFrame_jButtonPause_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonPause_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonPause_actionPerformed(e);
    }
}

class MainFrame_jButtonPreview_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonPreview_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonPreview_actionPerformed(e);
    }
}

class MainFrame_jButtonCancel_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonCancel_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonCancel_actionPerformed(e);
    }
}

class MainFrame_jButtonInit_actionAdapter
        implements ActionListener {

    MainFrame adaptee;

    MainFrame_jButtonInit_actionAdapter(MainFrame adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.jButtonInit_actionPerformed(e);
    }
}
