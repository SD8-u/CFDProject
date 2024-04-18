import sys
import os
import customtkinter
import textwrap
from mpi4py import MPI
from PIL import Image
from customtkinter import filedialog
from CTkMessagebox import CTkMessagebox

#Execute the simulation in parallel by initating the solver with MPI
def call_parallel_solver(refinement, timesteps, velocity, dt, viscosity, filename):
    comm = MPI.COMM_SELF.Spawn(
        sys.executable,
        args = [os.path.dirname(os.path.abspath(__file__)) + '/solver.py', 
                refinement, timesteps, velocity, dt, viscosity, filename],
        maxprocs=4
    )

    comm.barrier()

customtkinter.set_default_color_theme("dark-blue")

#Control window
class OptionFrame(customtkinter.CTkFrame):
    def __init__(self, master, image_update):
        super().__init__(master)
        self.master = master
        self.image_update = image_update
        title_font = customtkinter.CTkFont(family="<family name>", weight='bold')
        self.filepath = ""

        #Define UI widgets in option menu
        self.title = customtkinter.CTkLabel(self, text='CFD Simulation', font=title_font)
        self.button = customtkinter.CTkButton(self, text="Compute Solution", command=self.button_callbck)
        self.timestepsL = customtkinter.CTkLabel(self, text="Timesteps:")
        self.timesteps = customtkinter.CTkEntry(self, placeholder_text="Enter timesteps")
        self.viscosityL = customtkinter.CTkLabel(self, text="Viscosity:")
        self.viscosity = customtkinter.CTkEntry(self, placeholder_text="Enter viscosity")
        self.refinementL = customtkinter.CTkLabel(self, text="Refinement:")
        self.refinement = customtkinter.CTkEntry(self, placeholder_text="Enter refinement")
        self.dtL = customtkinter.CTkLabel(self, text="Time step size:")
        self.dt = customtkinter.CTkEntry(self, placeholder_text="Enter time step size")
        self.velocityL = customtkinter.CTkLabel(self, text="Velocity:")
        self.velocity = customtkinter.CTkEntry(self, placeholder_text="Enter velocity")
        self.file = customtkinter.CTkButton(self, text="Choose Geometry", command=self.file_explorer)
        self.filename = customtkinter.CTkLabel(self, text="File: ")

        #Position widgets in grid
        self.title.grid(row=0, column=0, padx=20, pady=5, sticky="w")
        self.button.grid(row=1, column=0, padx=20, pady=5, sticky="w")
        self.file.grid(row=2, column=0, padx=20, pady=5, sticky="w")
        self.filename.grid(row=3, column=0, padx=20, pady=5, sticky="w")
        self.timestepsL.grid(row=4, column=0, padx=20, pady=0, sticky="w")
        self.timesteps.grid(row=5, column=0, padx=20, pady=0, sticky="w")
        self.viscosityL.grid(row=6, column=0, padx=20, pady=5, sticky="w")
        self.viscosity.grid(row=7, column=0, padx=20, pady=0, sticky="w")
        self.refinementL.grid(row=8, column=0, padx=20, pady=5, sticky="w")
        self.refinement.grid(row=9, column=0, padx=20, pady=0, sticky="w")
        self.dtL.grid(row=10, column=0, padx=20, pady=5, sticky="w")
        self.dt.grid(row=11, column=0, padx=20, pady=0, sticky="w")
        self.velocityL.grid(row=12, column=0, padx=20, pady=5, sticky="w")
        self.velocity.grid(row=13, column=0, padx=20, pady=5, sticky="w")

        self.button.configure(state="disabled")

    #Execute simulation given parameters on button click
    def button_callbck(self):
        call_parallel_solver(self.refinement.get(), self.timesteps.get(), 
                           self.velocity.get(), self.dt.get(), self.viscosity.get(),
                           self.filepath)
        self.image_update()
    
    #Open file explorer to obtain geometry
    def file_explorer(self):
        filename = filedialog.askopenfilename(initialdir = os.getcwd() + "/geometry", 
                                              title = "Select a Geometry File", filetypes = 
                                              (("geo files", "*.geo*"), ("all files", "*.*")))
        self.filepath = filename
        filename = os.path.basename(filename)
        self.filename.configure(text = '\n'.join(textwrap.wrap('File: "'+filename+'"', 25)))
        self.button.configure(state="normal")

#Application window
class App(customtkinter.CTk):
    def __init__(self):
        #Initialise application UI
        super().__init__()
        customtkinter.set_appearance_mode("dark")
        #Display warning message prior to opening UI
        warning = CTkMessagebox(title="Warning", message="Solutions are merely an approximation - " + 
                                " do not expect complete accuracy.", icon="warning")
        warning.wait_window()
        self.title("CFD Solver")
        self.geometry("1200x680")
        self.image = customtkinter.CTkImage(dark_image=Image.open("stream_plot.png"), size=(940, 680))
        self.label = customtkinter.CTkLabel(master=self, image=self.image, text='')
        self.label.grid(row=0, column=1, padx=20, pady=20, sticky="w")
        self.option_frame = OptionFrame(self, self.image_update)
        self.option_frame.grid(row=0, column = 0, padx=20, pady=20, sticky="w")


    #Update image following simulation execution
    def image_update(self):
        self.image = customtkinter.CTkImage(dark_image=Image.open("stream_plot.png"), size=(940, 680))
        self.label.configure(image=self.image)

app = App()
app.mainloop()
