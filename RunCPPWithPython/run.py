import os


def run(f):
    cpp_f = f + '.cpp'
    os.system("echo Compiling " + cpp_f)
    os.system('g++ ' + cpp_f + ' -o ' + f)
    os.system("echo Running " + f)
    os.system("echo -------------------")
    os.system(f)


run('hello')
