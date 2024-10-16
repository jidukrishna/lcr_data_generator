import os
from tkinter import messagebox

from reportlab.lib import colors

import numpy as np
import matplotlib.pyplot as plt
from reportlab.platypus import SimpleDocTemplate, Table, Image, TableStyle, Paragraph, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter


def generate_graph(L,R,C,V,start_f,end_f,no_of_points,name_graph="graph.png"):



    def get_fr(L, C,R):
        """Function to compute resonant frequency."""
        f_resonant = 1 / (2 * np.pi * np.sqrt(L * C))
        f_r = (1 / (2 * np.pi)) * np.sqrt(1 / (L * C)) * np.sqrt(1 - (R / (2 * np.sqrt(L/C)))**2)

        return f_resonant


    def get_Imax(V, L, R, C):
        """Function to compute maximum current and resonant frequency."""
        f_resonant = get_fr(L, C,R)
        w = 2 * np.pi * f_resonant
        I_max = (V / np.sqrt(R**2 + (w * L - 1 / (w * C))**2)) * 10**3  # in mA
        return I_max, f_resonant


    # Define the frequency response function (amplitude of current I as a function of frequency)
    def current_amplitude(omega, V, R, L, C):
        return (V / np.sqrt(R**2 + (omega * L - 1 / (omega * C))**2)) * 10**3  # in mA


    # Get the maximum current and resonant frequency
    I_max, f_resonant = get_Imax(V, L, R, C)

    # Define the angular frequency range (in radians per second)
    frequencies = np.linspace(start_f, end_f, no_of_points)  # frequency in Hz
    angular_frequencies = 2 * np.pi * frequencies  # convert to angular frequency

    # Compute current amplitude over the range of frequencies
    I_omega = current_amplitude(angular_frequencies, V, R, L, C)

    # Calculate I_max * 0.707 (half-power point)
    i_707 = 0.707 * I_max

    # Function to find the frequencies where current equals I_half_power
    def find_half_power_frequencies(frequencies, I_omega, I_half_power):
        half_power_frequencies = []
        for i in range(len(I_omega) - 1):
            if (I_omega[i] - I_half_power) * (I_omega[i + 1] - I_half_power) < 0:
                # Linear interpolation to find exact frequency
                f_half = frequencies[i] + (frequencies[i + 1] - frequencies[i]) * (I_half_power - I_omega[i]) / (I_omega[i + 1] - I_omega[i])
                half_power_frequencies.append(f_half)
        return half_power_frequencies

    # Find the half-power frequencies
    half_power_frequencies = find_half_power_frequencies(frequencies, I_omega, i_707)

    # Plotting the frequency response
    fig, axe = plt.subplots()
    fig.set_figheight(6)
    fig.set_figwidth(10)
    axe.plot(frequencies, I_omega, label="Current Amplitude", color='b')





    def draw_lines():
        # Mark the resonant frequency and add it to the legend
        axe.axhline(I_max, color='orange', linestyle='--', label=f"Imax: {round(I_max, 2)} mA")

        # Draw horizontal line at I_max * 0.707 (half-power point) and add it to the legend
        axe.axhline(i_707, color='purple', linestyle='--', label=f"Imax*0.707: {round(i_707, 2)} mA")

        axe.axvline(x=f_resonant, color='r', linestyle='--', label=f"Resonant Frequency: {f_resonant:.2f} Hz")

        k = ("g", "black")
        c = 0
        g = "lower", "higher"
        for i, f_half in enumerate(half_power_frequencies):
            axe.axvline(x=f_half, color=k[c], linestyle='--', label=f"{g[c]} Frequency {i + 1}: {f_half:.2f} Hz")
            c += 1




    def scatter_plots():
        axe.scatter(f_resonant, I_max, color='r', label=f"{f_resonant:.2f} Hz {I_max:.2f} mA")
        k = ("g", "black")
        c = 0
        for i, f_half in enumerate(half_power_frequencies):
            axe.scatter(f_half, i_707, color=k[c], label=f"{f_half:.2f} Hz {i_707:.2f} mA")
            c += 1

    draw_lines()
    scatter_plots()





    axe.set_title("Frequency Response of LCR Circuit")
    axe.set_xlabel("Frequency [Hz]")
    axe.set_ylabel("Current Amplitude [mA]")
    axe.grid(True)
    axe.legend(loc='upper right')  # Adjust legend location

    plt.savefig(name_graph)

    return (half_power_frequencies[0].item(),half_power_frequencies[1].item(),
            f_resonant.item(),I_max.item(), i_707.item(),frequencies,I_omega)





def generate_pdf(L, R, C, V, start_f, end_f, no_of_points,pdf_file=f"electrical damping lcr.pdf",image_name_graph="graph.png"):
    k = generate_graph(L, R, C, V, start_f, end_f, no_of_points, image_name_graph)
    pdf = SimpleDocTemplate(pdf_file, pagesize=letter)

    styles = getSampleStyleSheet()


    # Data for the table
    heading = ['frequency hz', 'Current (mA)']

    freq, currents = k[-2:]
    data = [heading]
    max_current = 0
    c = 0
    j = 0
    for i in range(len(freq)):
        c += 1
        data.append([freq[i], currents[i]])
        if currents[i] > max_current:
            max_current = currents[i]
            j = c


    table = Table(data)

    # Add some style to the table
    style = TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.grey),  # Header background
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),  # Header text color
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),  # Align all cells center
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),  # Header font
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),  # Padding for header
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),  # Background color for body
        ('GRID', (0, 0), (-1, -1), 1, colors.black),  # Gridlines
    ])

    table.setStyle(style)
    table.setStyle([('BACKGROUND', (0, j), (-1, j), colors.violet)])


    #image
    image_name = Image(image_name_graph)
    res = 0.7
    image_name.drawWidth = 1000 * res
    image_name.drawHeight = 600 * res



    #paragraphs

    main_heading = Paragraph(f"electrical damping for an LCR circuit", styles["Heading1"])


    para_list = []

    b = [

        f"""L = {L} H, R= {R} ohm C = {C} F , V = {V} V, """,
        f"""I max = {round(k[3], 3)} mA""",
        f"""I max * 0.707 = {round(k[4], 3)} mA""",

        f"""
    Resonant frequency = {round(k[2], 3)} Hz\n

        """, f"""
    Resonant frequency = {round(k[2], 3)} Hz\n

        """, f'''

    Band width= fh -fl = {round(k[1], 3)} - {round(k[0], 3)} hz

    ''', f'''

    Band width = {round(k[1] - k[0], 3)} Hz\n
    ''',
        f'''

    Quality factor = fr/band width = {round(k[2], 3)}/{round(k[1] - k[0], 3)}
    ''', f'''

    Quality factor =  {round(k[2] / (k[1] - k[0]), 3)}
    ''',
    ]

    for i in b:
        para_list.append(Paragraph(i, styles['Normal']))

    pdf.build([main_heading, table, image_name, *para_list])









import customtkinter as c

def main():
    # LCR circuit parameters
    L = 30 * 10 ** -3  # inductance in Henry
    R = 160  # resistance in Ohms
    C = 0.1 * 10 ** -6  # capacitance in Farads
    V = 5.0  # voltage source amplitude in Volts
    start_f = 1  # starting frequency in Hz
    end_f = 10000  # ending frequency in Hz
    no_of_points = 100  # number of points for frequency sweep
    name_of_graph_pic = "graph.png"
    pdf_file = f"electrical damping lcr.pdf"
    #generate_pdf(L, R, C, V, start_f, end_f, no_of_points, pdf_file,name_of_graph_pic)
    root=c.CTk()
    root.title("LRC CIRCUIT DATA GENERATOR")
    root.geometry("600x450")
    c.CTkLabel(root,text="LRC CIRCUIT DATA GENERATOR").grid(column=0, row=0)
    c.CTkLabel(root,text="L (H) :").grid(column=0, row=1)
    c.CTkLabel(root,text="R (ohm) :").grid(column=0, row=2)
    c.CTkLabel(root,text="C (F) :").grid(column=0, row=3)
    c.CTkLabel(root,text="V (V) :").grid(column=0, row=4)
    c.CTkLabel(root,text="start frequency (hz) :").grid(column=0, row=5)
    c.CTkLabel(root,text="end frequency (hz) :").grid(column=0, row=6)
    c.CTkLabel(root,text="no of plots for graph :").grid(column=0, row=7)
    c.CTkLabel(root,text="name of pdf :").grid(column=0, row=8)
    c.CTkLabel(root,text="name of graph :").grid(column=0, row=9)
    c.CTkLabel(root,text="1. fill all the columns").grid(column=0, row=11)
    c.CTkLabel(root,text="2. for multiplication use this format 3*10**2\nwhich is same as 3x10^2").grid(column=0, row=12,rowspan=2)
    c.CTkLabel(root,text="3. use .pdf and .png while saving").grid(column=0, row=14)
    c.CTkLabel(root,text="made by jidu krishna , cs-c 28\'").grid(column=1, row=15)


    L_entry = c.CTkEntry(root,width=250)
    R_entry = c.CTkEntry(root,width=250)
    C_entry = c.CTkEntry(root,width=250)
    V_entry = c.CTkEntry(root,width=250)
    start_f_entry = c.CTkEntry(root,width=250)
    end_f_entry = c.CTkEntry(root,width=250)
    no_of_points_entry = c.CTkEntry(root,width=250)
    pdf_file_entry = c.CTkEntry(root,width=250)
    graph_name_entry = c.CTkEntry(root,width=250)

    L_entry.grid(column=1, row=1)
    R_entry.grid(column=1, row=2)
    C_entry.grid(column=1, row=3)
    V_entry.grid(column=1, row=4)
    start_f_entry.grid(column=1, row=5)
    end_f_entry.grid(column=1, row=6)
    no_of_points_entry.grid(column=1, row=7)
    pdf_file_entry.grid(column=1, row=8)
    graph_name_entry.grid(column=1, row=9)



    def make_it():
        V=float(eval(V_entry.get()))
        R=float(eval(R_entry.get()))
        C=float(eval(C_entry.get()))
        L=float(eval(L_entry.get()))
        start_f=int(start_f_entry.get())
        end_f=int(end_f_entry.get())
        no_of_points=int(no_of_points_entry.get())
        pdf_file=pdf_file_entry.get()
        name_of_graph_pic=graph_name_entry.get()
        generate_pdf(L, R, C, V, start_f, end_f, no_of_points, pdf_file,name_of_graph_pic)
        messagebox.showinfo("Successful", f"Data generation successful and stored in {os.getcwd()}/{pdf_file}")


    submit=c.CTkButton(root,text="Submit",command=lambda:make_it())
    submit.grid(column=0, row=10,columnspan=2,padx=5,pady=5)













    root.mainloop()


if __name__ == "__main__":
    main()
