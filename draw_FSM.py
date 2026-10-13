import graphviz

dot = graphviz.Digraph('MooreMachine', engine='circo', filename='3D_Printer_FSM_Circular_Tweaked')

dot.attr(overlap='false', splines='true')
dot.attr('node', shape='Mrecord', fontname='Times-Roman', fontsize='12', margin='0.2')
dot.attr('edge', fontname='Times-Roman', fontsize='11', len='4.0')

dot.node('S0', '{ S0: Printing | \u0393: Extrude_Normal }')
dot.node('S1', '{ S1: Stop_Save | \u0393: Save_XY }')
dot.node('S2', '{ S2: Moving | \u0393: Move_To_0_0 }')
dot.node('S3', '{ S3: Wiping | \u0393: Extrude_Wipe }')
dot.node('S4', '{ S4: Cooling | \u0393: Heat_OFF }')
dot.node('S5', '{ S5: Loading | \u0393: Feed_Filament }')
dot.node('S6', '{ S6: Heating | \u0393: Heat_ON }')
dot.node('S7', '{ S7: Returning | \u0393: Move_To_Saved_XY }')

dot.node('start', 'Start', shape='plaintext', fontname='Times-Roman', fontsize='12')
dot.edge('start', 'S0', arrowsize='1.5', penwidth='2.0', minlen='1.5')

dot.edge('S0', 'S1', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>C (New_Colour)<BR/>OR F (Filament_Finished)</TD></TR></TABLE>>')
dot.edge('S1', 'S2', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>&nbsp;&nbsp;&nbsp;Position_Saved<BR/>&nbsp;</TD></TR></TABLE>>')
dot.edge('S2', 'S3', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>&nbsp;<BR/>&nbsp;&nbsp;&nbsp;Parked_At_Origin</TD></TR></TABLE>>')
dot.edge('S3', 'S4', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>&nbsp;<BR/>Excess_Removed</TD></TR></TABLE>>')
dot.edge('S4', 'S5', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>Temp_Low</TD></TR></TABLE>>')
dot.edge('S5', 'S6', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>&nbsp;&nbsp;&nbsp;Filament_Loaded<BR/>&nbsp;</TD></TR></TABLE>>')
dot.edge('S6', 'S7', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>Temp_High&nbsp;&nbsp;&nbsp;<BR/>&nbsp;</TD></TR></TABLE>>')
dot.edge('S7', 'S0', label='<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="2" BGCOLOR="white"><TR><TD>&nbsp;<BR/>Returned_To_Saved</TD></TR></TABLE>>')

dot.render(view=True, format='pdf', cleanup=True)
print("Tweaked diagram generated successfully as a PDF.")