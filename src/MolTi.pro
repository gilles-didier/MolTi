#-------------------------------------------------
#
# Project created by QtCreator 2015-04-03T15:59:33
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets
QT       += webkitwidgets

TARGET = MolTi
TEMPLATE = app
ICON = MolTi.icns


SOURCES += main.cpp\
		moltiwindow.cpp \
	Annotation.c \
	EdgesComposition.c \
	Louvain.c \
	Partition.c \
	StatAnnotationPartition.c \
	Utils.c \
	Graph.c \
	qtablegraphmodel.cpp \
	qhelpdialog.cpp \
	qpartitionmodel.cpp \
	qpartitionvaluemodel.cpp \
	qitemclassdelegate.cpp \
	qitemdoubledelegate.cpp \
	qpartitionitemdelegate.cpp \
	qitemintdelegate.cpp \
	qitembigstringdelegate.cpp

HEADERS  += moltiwindow.h \
	Annotation.h \
	EdgesComposition.h \
	Louvain.h \
	Partition.h \
	StatAnnotationPartition.h \
	Utils.h \
	Graph.h \
	qtablegraphmodel.h \
	qhelpdialog.h \
	qpartitionmodel.h \
	qpartitionvaluemodel.h \
	qitemclassdelegate.h \
	qitemdoubledelegate.h \
	qpartitionitemdelegate.h \
	qitemintdelegate.h \
	qitembigstringdelegate.h

FORMS    += moltiwindow.ui \
	qhelpdialog.ui

DISTFILES += \
	MolTi.pro.user

unix|win32: LIBS += -lgsl

unix|win32: LIBS += -lgslcblas

RESOURCES += \
	data/resources.qrc
