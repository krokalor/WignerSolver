# WignerSolver
Program został napisany w języku C++ i oparty jest na bibliotece numerycznej armadillo.
Część obliczeń oparta jest na macierzach rzadkich, więc koniecznym jest też zainstalowanie bilbioteki SuperLU.
#
Pliki zawierające deklaracje klas oraz funkcji znajdują się w folderze 'src'.
Program obsluguje się z poprzez plik 'main.cpp'.
W celu skompilowania programu należy wpisać komendę 'make' będąc w folderze WignerSolver.
Żeby wyłączyć zrównoleglanie należy zakomentować w pliku makefile fragment '-fopenmp'.
