#!/usr/bin/env python3

import shutil, os, argparse, sys, stat, subprocess, shlex
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyUtils"))
from utils import Utils

class genHelper:
    # the C++ standard the projects are compiled with; the compiler check below makes sure the chosen compiler
    # supports it, rather than relying on a versioned compiler name (e.g. g++-10) to imply a minimum version
    cxxStandard = "c++17"

    @staticmethod
    def generateCompfileFull(outFileName, externalDirLoc, cc, cxx, outName, installDirName, installDirLoc, neededLibs,ldFlags = "", cxxFlags = "", private =False):
        availableLibs = ["CPPITERTOOLS","CPPPROGUTILS","ZI_LIB","BOOST",
                         "R","BAMTOOLS","CPPCMS","MATHGL","ARMADILLO",
                         "MLPACK","LIBLINEAR","CURL","GTKMM", "NJHSEQ",
                          "NJHCPP", "SEEKDEEP",
                         "NJHSEQDEV", "SEEKDEEPDEV", "CATCH", "CATCH2", "JSONCPP",
                          "TWOBIT", "SEQSERVER","NJHRINSIDE", "PSTREAMS",
                           "MONGOC", "MONGOCXX", "SHAREDMUTEX",
                           "MAGIC", "HTS", "RESTBED", "LIBPCA", "BOOST_FILESYSTEM","BOOST_MATH",
                           "PATHWEAVER", "ELUCIDATOR", "ELUCIDATORLAB", "MIPWRANGLER", "EIGEN",
                            "ZLIB", "ZLIB-NG", "UNQLITE", "NLOHMANN_JSON"]
        neededLibraries = {}
        for lib in neededLibs:
            if ":" in lib:
                libSplit = lib.split(":")
                neededLibraries[libSplit[0].upper()] = ":".join(libSplit[1:])
            else:
                neededLibraries[lib.upper()] = ""

        with open(outFileName, "w") as f:
            f.write("CC = {CC}\n".format(CC = cc))
            f.write("CXX = {CXX}\n".format(CXX = cxx))
            f.write("CXXOUTNAME = {NAME_OF_PROGRAM}\n".format(NAME_OF_PROGRAM = outName))
            f.write("CXXFLAGS = -std=" + genHelper.cxxStandard + "\n")
            #f.write("CXXFLAGS += -Wall -ftemplate-depth=1024 -Werror=uninitialized -Werror=return-type -Wno-missing-braces\n")
            f.write("CXXFLAGS += -Wall -ftemplate-depth=1024 -Wno-deprecated-declarations -Werror=uninitialized -Werror=return-type -Wno-missing-braces\n")
            if "" != cxxFlags:
                if cxxFlags.startswith("\\"):
                    cxxFlags = cxxFlags[1:]
                f.write("CXXFLAGS += " + cxxFlags +"\n")
            if "" != ldFlags:
                f.write("LD_FLAGS = ")
                if not ldFlags.startswith("-"):
                    f.write("-")
                f.write("{ld_flags}\n".format(ld_flags = " ".join(ldFlags.split(","))))
            f.write("CXXOPT += -O2 -funroll-loops -DNDEBUG  \n")
            f.write("\n")
            if private:
                f.write("PRIVATE = TRUE\n")
            f.write("#debug\n")
            f.write("CXXDEBUG = -g -gstabs+ \n")
            f.write("INSTALL_DIR={INSTALL_LOCATION}\n".format(INSTALL_LOCATION = os.path.join(installDirLoc,installDirName)))
            f.write("EXT_PATH=$(realpath {EXTERNAL})\n".format(EXTERNAL = externalDirLoc))
            f.write("\n")
            for lib in availableLibs:
                if lib in neededLibraries:
                    if neededLibraries[lib] == "":
                        f.write("USE_{LIB} = 1\n".format(LIB = lib))
                    else:
                        f.write("USE_{LIB} = 1#{BRANCH}\n".format(LIB = lib, BRANCH = neededLibraries[lib]))
                #else:
                #    f.write("USE_{LIB} = 0\n".format(LIB = lib))



    @staticmethod
    def compilerVersion(compiler):
        """First line of `compiler --version`, or "" if it can't be run."""
        try:
            out = subprocess.run(shlex.split(compiler) + ["--version"], stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, universal_newlines=True, timeout=60).stdout
            return out.strip().splitlines()[0] if out.strip() else ""
        except (OSError, subprocess.SubprocessError):
            return ""

    @staticmethod
    def checkCxxStandardSupport(cxx, standard = None):
        """Test compile a small program using the C++17 features these projects rely on.
        Returns (True, "") if cxx can compile it, otherwise (False, the compiler's output)."""
        if standard is None:
            standard = genHelper.cxxStandard
        src = """
#if __cplusplus < 201703L
#error "the compiler does not support C++17"
#endif
#include <any>
#include <optional>
#include <string_view>
#include <variant>
#include <utility>
int main() {
  auto [a, b] = std::pair<int, int>{1, 2};
  if constexpr (sizeof(int) > 0) { a += 1; }
  std::optional<std::string_view> o{"ok"};
  std::variant<int, double> v{2.0};
  std::any x = a + b;
  return (o && std::holds_alternative<double>(v) && x.has_value()) ? 0 : 1;
}
"""
        try:
            proc = subprocess.run(shlex.split(cxx) + ["-std=" + standard, "-fsyntax-only", "-x", "c++", "-"],
                                  input=src, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  universal_newlines=True, timeout=120)
        except (OSError, subprocess.SubprocessError) as e:
            return False, str(e)
        return proc.returncode == 0, proc.stdout

    @staticmethod
    def firstCompilerErrors(output, n = 5):
        """The first few error lines of compiler output, or its first lines if none say error."""
        lines = output.strip().splitlines()
        errors = [l for l in lines if "error" in l]
        return "\n".join((errors or lines)[:n])

    @staticmethod
    def requireCxxStandardSupport(cxx, standard = None):
        """Exit with an explanation if cxx can't compile the C++ standard these projects need."""
        if standard is None:
            standard = genHelper.cxxStandard
        ok, output = genHelper.checkCxxStandardSupport(cxx, standard)
        if ok:
            return
        version = genHelper.compilerVersion(cxx)
        msg = "Error: the C++ compiler " + cxx + (" (" + version + ")" if version else "") + \
              " can't compile " + standard + ", which is needed.\n"
        if not version:
            msg += "It could not be run, check that it is installed and on your PATH.\n"
        else:
            msg += "Compiler output:\n" + genHelper.firstCompilerErrors(output) + "\n"
        msg += "Give a newer compiler with -CC and -CXX or the CC and CXX environment variables, " + \
               "e.g. CC=gcc-13 CXX=g++-13 (gcc 7 or newer, or clang 5 or newer, is needed for C++17).\n"
        sys.stderr.write(msg)
        sys.exit(1)

    @staticmethod
    def determineCC(args, defaultCC = "gcc"):
        if Utils.isMac():
            defaultCC = "clang"
        if not args.CC:
            eCC = os.getenv("CC")
            if(eCC):
                defaultCC =  eCC
        else:
            defaultCC =  Utils.getStrFromStrOrList(args.CC)
        return defaultCC

    @staticmethod
    def determineCXX(args, defaultCXX = "g++"):
        if Utils.isMac():
            defaultCXX = "clang++"
        if not args.CXX:
            eCXX = os.getenv("CXX")
            if  eCXX:
                defaultCXX = eCXX
        else:
            defaultCXX = Utils.getStrFromStrOrList(args.CXX)
        return defaultCXX

    @staticmethod
    def parseNjhConfigureArgs():
        parser = argparse.ArgumentParser()
        parser.add_argument('-prefix', type=str, nargs=1)
        parser.add_argument('-externalLibDir', type=str, nargs=1)
        parser.add_argument('-CC', type=str, nargs=1)
        parser.add_argument('-CXX', type=str, nargs=1)
        args = parser.parse_args()
        # check here, in configure.py's own process, so a compiler that's too old gives a plain message
        # rather than one wrapped in the tracebacks of the scripts configure.py runs
        genHelper.requireCxxStandardSupport(genHelper.determineCXX(args))
        return args

    @staticmethod
    def mkConfigCmd(name,libs, argv, ldflags="", cxxFlags="", private = False):
        if libs == "":
            cmd = os.path.join(os.path.dirname(os.path.dirname(__file__)), "setUpScripts/njhConfigure.py") + " -name {name} ".format(name=name)
        else:
            cmd = os.path.join(os.path.dirname(os.path.dirname(__file__)), "setUpScripts/njhConfigure.py") + " -name {name} -libs {libs}".format(name=name, libs=libs)
        if "" != ldflags:
            if ldflags.startswith("-"):
                ldflags = ldflags[1:]
            cmd += " -ldFlags " + ldflags
        if "" != cxxFlags:
            addingFlags = " -cxxFlags \""
            if cxxFlags.startswith("-"):
                addingFlags += "\\"
            cmd += addingFlags + cxxFlags + "\""
        if len(sys.argv) > 1:
            cmd += " " + " ".join(sys.argv[1:])
        if private:
            cmd += " -private "
        return cmd


    @staticmethod
    def mkConfigFileStr(name, libs):
        ret = """#!/usr/bin/env python2
import shutil, os, argparse, sys, stat
sys.path.append("scripts/pyUtils")
sys.path.append("scripts/setUpScripts")
from utils import Utils
from genFuncs import genHelper

def main():
    name = "{name}"
    libs = "{libs}"
    args = genHelper.parseNjhConfigureArgs()
    cmd = genHelper.mkConfigCmd(name, libs, sys.argv)
    Utils.run(cmd)

main()
""".format(name = name, libs = libs)
        return ret
