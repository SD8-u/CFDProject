import customtkinter

customtkinter.set_default_color_theme("dark-blue")

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()
        customtkinter.set_appearance_mode("dark")
        self.title("Blood Flow CFD")
        self.geometry("1000x800")
        self.grid_columnconfigure((0, 1), weight=1)
        self.button = customtkinter.CTkButton(self, text="compute solution", command=self.button_callbck)
        self.button.pack(padx=20, pady=20)

    def button_callbck(self):
        print("Solve system")

app = App()
app.mainloop()