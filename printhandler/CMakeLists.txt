cmake_minimum_required(VERSION 2.8)

project(PrintHandler)

set(srcs
	SRC/TC_Common/RegExp.cxx
	SRC/TC_Common/TC_Bit.cxx
	SRC/TC_Common/TC_MinGrid.cxx
	SRC/TC_Common/TC_Name.cxx
	SRC/TC_Common/TC_NameLength.cxx
	SRC/TC_Common/TC_NameType.cxx
	SRC/TC_Common/TC_SideIndex.cxx
	SRC/TC_Common/TC_SideName.cxx
	SRC/TC_Common/TC_StringUtils.cxx
	SRC/TIO_InputOutputHandlers/TIO_CustomOutput.cxx
	SRC/TIO_InputOutputHandlers/TIO_FileHandler.cxx
	SRC/TIO_InputOutputHandlers/TIO_FileOutput.cxx
	SRC/TIO_InputOutputHandlers/TIO_PrintHandler.cxx
	SRC/TIO_InputOutputHandlers/TIO_PrintHandlerExtern.cxx
	SRC/TIO_InputOutputHandlers/TIO_SkinHandler.cxx
	)

add_library(PrintHandler STATIC
	${srcs}
)

target_include_directories(PrintHandler PUBLIC
	../pcre/SRC
	SRC/TC_Common
	SRC/TIO_InputOutputHandlers
)

target_compile_options(PrintHandler PRIVATE
	-ansi -pedantic
	-Wpointer-arith -Wcast-qual -Wshadow -Wcast-align -Wno-write-strings
	-D__USE_FIXED_PROTOTYPES__ -D_POSIX_SOURCE -DLINUX24_64 -DLINUX_X86_64 -D EZXML_NOMMAP -D_POSIX_C_SOURCE)
