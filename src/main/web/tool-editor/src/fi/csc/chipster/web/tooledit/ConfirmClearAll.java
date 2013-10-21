package fi.csc.chipster.web.tooledit;

import com.vaadin.ui.Alignment;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window;

/**
 * Confirmation window for clearing all data
 * @author Gintare Pacauskaite
 *
 */
public class ConfirmClearAll extends Window{

	private static final long serialVersionUID = -6901652612577754156L;
	private Button ok = new Button("Ok");
	private Button cancel = new Button("Cancel");
	private ToolEditorUI root;

	public ConfirmClearAll(ToolEditorUI root) {
		super("Do you really want to clear everyting?");
		this.setModal(true);
		this.setHeight("200px");
		this.root = root;
		initButtonsListener();
		VerticalLayout vLayout = new VerticalLayout();
		vLayout.setSpacing(true);
		vLayout.setMargin(true);
		vLayout.addComponent(new Label("This will clear all inputs, outputs, and parameters as well as the the text area. Do you really want to clear all?"));
		HorizontalLayout hLayout = new HorizontalLayout();
		hLayout.addComponent(ok);
		hLayout.addComponent(cancel);
		hLayout.setSpacing(true);
		vLayout.addComponent(hLayout);
		vLayout.setComponentAlignment(hLayout, Alignment.BOTTOM_RIGHT);
		this.setContent(vLayout);
	}
	
	private void initButtonsListener() {
		ok.addClickListener(new ClickListener() {
			private static final long serialVersionUID = 6970174755245561409L;

			@Override
			public void buttonClick(ClickEvent event) {
				root.getToolEditor().removeAllComponents();
				root.getTreeToolEditor().removeAllChildren();
				root.getTreeToolEditor().updateToolTitle("");
				root.getTextEditor().clearAllText();
				root.removeWindow(ConfirmClearAll.this);
			}
		});
		
		cancel.addClickListener(new ClickListener() {
			private static final long serialVersionUID = 5621280789210184012L;

			@Override
			public void buttonClick(ClickEvent event) {
				root.removeWindow(ConfirmClearAll.this);
			}
		});
	}

}
