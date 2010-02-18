package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.FocusEvent;

import javax.swing.JComponent;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

public class StringInputComponent extends ParameterInputComponent implements CaretListener, DocumentListener { 

	private StringParameter param;
	private JTextField field;

	protected StringInputComponent(StringParameter parameter, ParameterPanel parameterPanel) {		
		super(parameterPanel);
		this.param = parameter;
		this.field = new JTextField();
		field.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);
		field.setText("" + parameter.getValue());
		field.addCaretListener(this);
		field.addFocusListener(this);
		field.getDocument().addDocumentListener(this);
		this.add(field, BorderLayout.CENTER);
	}

	
	@Override
	public Parameter getParameter() {
		return param;
	}

	@Override
	public boolean inputIsValid() {
		return true;
	}

	public void caretUpdate(CaretEvent e) {
		getParentPanel().setMessage(param.getDescription(), Color.black);
	}

	private void updateParameter() {
		param.setValue(field.getText());
	}
		
	public void changedUpdate(DocumentEvent e) {
		updateParameter();
	}

	public void insertUpdate(DocumentEvent e) {
		updateParameter();
	}

	public void removeUpdate(DocumentEvent e) {
		updateParameter();
	}


	@Override
	public JComponent getParameterComponent() {
		return field;
	}

	public void focusGained(FocusEvent e) {
		getParentPanel().setMessage(param.getDescription(), Color.black);
	}

	
}
