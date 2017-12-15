% Klasa dokumentow do skladu prac inżynierskich, licencjackich i magisterskich na wydziale MiNI PW.
%------------------------------------------------------------------------------%
\NeedsTeXFormat{LaTeX2e}[2016/10/25]
\ProvidesClass{minipw}[2016/10/25 Szablon pracy dyplomowej MiNI PW]
%------------------------------------------------------------------------------%
\RequirePackage{ifthen}
\newboolean{lady}
\setboolean{lady}{false}

\newif\ifmyclass@en
\DeclareOption{en}{\myclass@entrue}
\DeclareOption{pl}{\myclass@enfalse}
\ExecuteOptions{pl}
\ProcessOptions\relax

\DeclareOption{lady}{\setboolean{lady}{true}}
\DeclareOption{11pt}{\OptionNotUsed}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax

% ---------------------------------------------------------

\LoadClass[a4paper,11pt,twoside]{report}
%\LoadClass[a4paper,11pt,twoside,openright]{report}
% to drugie gdyby ktoś chciał, żeby rozdziały się zaczynały zawsze na nieparzystej stronie

%------------------------------------------------------------------------------
\RequirePackage{graphicx}
\graphicspath{{graphics/}}

%\RequirePackage{caption} % jeśli używamy xelatexa, jest potrzebny, żeby robiły się \caption pod obrazkami i tabelami

%------------------------------------------------------------------------------%
\RequirePackage{geometry}
\geometry{
  inner=20mm,
  outer=20mm,
  bindingoffset=10mm,
  top=25mm,
  bottom=25mm,  % 20mm
  %includeheadfoot
}

\linespread{1.5}
%interlinia, zalecana przez Rektora 1.15, ale często kiedy w tekście jest dużo symboli typu suma lub całka, to one się po prostu nie mieszczą w linii; dlatego domyślnie 1.5 (ale więcej to już przesada, standard to 1.2 - 1.5)

% localsize pozwala zmieniać rozmiar czcionki w wybranym  środowisku
\newenvironment{localsize}[1]
{%
  %\clearpage
  \let\orignewcommand\newcommand
  \let\newcommand\renewcommand
  \makeatletter
  \input{bk#1.clo}%
  \makeatother
  \let\newcommand\orignewcommand
}{}

%------------------------------------------------------------------------------% Żywa pagina = napis na górze (parzysta - obecna sekcja, nieparzysta - obecny rozdział) + napis na dole (numer strony)

\RequirePackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
% numery stron: lewa do lewego, prawa do prawego 
\fancyfoot[LE,RO]{\thepage} 
% prawa pagina: zawartość \rightmark do lewego, wewnętrznego (marginesu) 
\fancyhead[LO]{\sc \nouppercase{\rightmark}}
% lewa pagina: zawartość \leftmark do prawego, wewnętrznego (marginesu) 
\fancyhead[RE]{\sc \leftmark}

\renewcommand{\chaptermark}[1]{
\markboth{\thechapter.\ #1}{}}

% kreski oddzielające paginy (górną i dolną):
\renewcommand{\headrulewidth}{0 pt} % 0 - nie ma, 0.5 - jest linia


\fancypagestyle{plain}{% to definiuje wygląd pierwszej strony nowego rozdziału - obecnie tylko numeracja
  \fancyhf{}%
  \fancyfoot[LE,RO]{\thepage}%
  
  \renewcommand{\headrulewidth}{0pt}% Line at the header invisible
  \renewcommand{\footrulewidth}{0.0pt}
}


%------------------------------------------------------------------------------
%Chaptery

\RequirePackage{titlesec}
\titleformat{\chapter}%[display]
  {\normalfont\Large \bfseries}
  {\thechapter.}{1ex}{\Large}

\titleformat{\section}
  {\normalfont\large\bfseries}
  {\thesection.}{1ex}{}
\titlespacing{\section}{0pt}{30pt}{20pt} 
%\titlespacing{\co}{akapit}{ile przed}{ile po} 
  
  
\titleformat{\subsection}
  {\normalfont \bfseries}
  {\thesubsection.}{1ex}{}
%------------------------------------------------------------------------------% Komenda \cleardoublepage powoduje, że to, co zostanie napisane po niej pojawi się na następnej nieparzystej stronie; czyli jeśli wywołujemy ją na parzystej, to przeskakuje do następnej, a jak na nieprzystej, robi jedną pustą stronę przerwy.

\def\cleardoublepage{\clearpage\if@twoside
\ifodd\c@page\else\hbox{}\thispagestyle{empty}\newpage
\if@twocolumn\hbox{}\newpage\fi\fi\fi}

%------------------------------------------------------------------------------%
%\RequirePackage[parfill]{parskip}[2001/04/09] % ten pakiet, jeśli jest włączony, umożliwia ustawienie odstępu po akapicie - jak w Wordzie

%------------------------------------------------------------------------------%
\clubpenalty=10000
\widowpenalty=10000
\brokenpenalty=10000
\hyphenation{Syl-ves-tra}
\hyphenation{Syl-ves-ter-a}
%------------------------------------------------------------------------------%
\RequirePackage{amsmath}
\RequirePackage{amsthm}
\RequirePackage{amsfonts}
\RequirePackage{mathrsfs} 
\RequirePackage{amssymb}


%------------------------------------------------------------------------------%
\RequirePackage[utf8]{inputenc}
\RequirePackage[MeX]{polski}
\RequirePackage[T1]{fontenc}
\RequirePackage[english, polish]{babel}

%------------------------------------------------------------------------------%
\RequirePackage[unicode]{hyperref}
\hypersetup{
  bookmarksopen,
  bookmarksopenlevel=1,
  pdfborder={0 0 0},
  pdfdisplaydoctitle,
  pdfpagelayout={SinglePage},
  pdfstartview={FitH}
}

%------------------------------------------------------------------------------%Spis treści

% wyglad spisu tresci
\renewcommand{\tableofcontents}{
\if@openright\cleardoublepage\else\clearpage\fi\thispagestyle{empty}
\pdfbookmark[0]{\contentsname}{tableofcontents}
\@makeschapterhead{\contentsname}\@afterheading
\@mkboth{\MakeUppercase\contentsname}{\MakeUppercase\contentsname}
\@starttoc{toc}}


% kropki dla chapterów
\RequirePackage{etoolbox} 
\patchcmd{\l@chapter}
  {\hfil}
  {\leaders\hbox{\normalfont$\m@th\mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}\hfill}
  {}{}


% Uwzględnia nienumerowane rozdzialy w spisie tresci
\def\@schapter#1{
\phantomsection\addcontentsline{toc}{chapter}{#1}
\markboth{#1}{#1}\@makeschapterhead{#1}\@afterheading}

% Stawia kropki po numerach rozdziałów, sekcji, podsekcji w spisie treści:

\RequirePackage{titletoc} 
\titlecontents{chapter}% <section-type>
  [0pt]% <left>
  {}% <above-code>
  {\bfseries \thecontentslabel.\quad}% <numbered-entry-format>
  {\bfseries}% <numberless-entry-format>
  {\bfseries\leaders\hbox{\normalfont$\m@th\mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}\hfill\contentspage}% <filler-page-format>

\titlecontents{section}
  [1em]
  {}
  {\thecontentslabel.\quad}
  {}
  {\leaders\hbox{\normalfont$\m@th\mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}\hfill\contentspage}

\titlecontents{subsection}
  [2em]
  {}
  {\thecontentslabel.\quad}
  {}
  {\leaders\hbox{\normalfont$\m@th\mkern \@dotsep mu\hbox{.}\mkern \@dotsep mu$}\hfill\contentspage}

%-------------------------------------------------------

% Ustawienia numeracji tabel i obrazkow
\renewcommand*{\thetable}{\arabic{chapter}.\arabic{table}}
\renewcommand*{\thefigure}{\arabic{chapter}.\arabic{figure}}
%\let\c@table\c@figure % jeśli włączone, numeruje tabele i obrazki razem

%------------------------------------------------------------------------------%Definicje

\def\title#1{\gdef\@title{#1}\hypersetup{pdftitle={#1}}}
\def\@title{\@latex@error{No title given}\@ehc}
\def\titleaux#1{\gdef\@titleaux{#1}}
\def\@titleaux{\@latex@error{No polish title given}\@ehc}
\def\author#1{\gdef\@author{#1}\hypersetup{pdfauthor={#1}}}
\def\@author{\@latex@error{No author given}\@ehc}
\def\supervisor#1{\gdef\@supervisor{#1}}
\def\@supervisor{\@latex@error{No supervisor given}\@ehc}

\def\konsultacje#1{\gdef\@konsultacje{#1}}
\def\@konsultacje{\@latex@error{Nie ma podanych konsultacji}\@ehc}

\def\type#1{\gdef\@type{#1}}
\def\@type{\@latex@error{No type given}\@ehc}
\def\discipline#1{\gdef\@discipline{#1}}
\def\@discipline{\@latex@error{No discipline given}\@ehc}
\def\specjal#1{\gdef\@specjal{#1}}
\def\@specjal{\@latex@error{Nie podano specjalnosci}\@ehc}
\def\date#1{\gdef\@date{#1}}
\def\@date{\@latex@error{No date given}\@ehc}
\def\album#1{\gdef\@album{#1}}
\def\@album{\@latex@error{No album number given}\@ehc}
\def\keywords#1{\gdef\@keywords{#1}}
\def\@keywords{\@latex@error{No keywords given}\@ehc}
\def\klucze#1{\gdef\@klucze{#1}}
\def\@klucze{\@latex@error{Nie podano slow kluczowych}\@ehc}

%-----------------------------------------------------------------------------

%------------------------------------------------------------------------------%
\newtheoremstyle{thm}{3ex}{3ex}{}{}{\bfseries}{.}{.5em}{#1 #2\thmnote{ (#3)}}
\newtheoremstyle{rem}{3ex}{3ex}{}{}{\bfseries}{.}{.5em}{#1 #2\thmnote{ (#3)}}
\theoremstyle{thm}



\ifmyclass@en

\theoremstyle{thm}
\newtheorem{theorem}{Theorem}[chapter]
\newtheorem{definition}[theorem]{Definition}
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{remark}[theorem]{Remark}
\newtheorem{example}[theorem]{Example}

\else
\theoremstyle{thm}
\newtheorem{theorem}{Twierdzenie}[chapter]
\newtheorem{definition}[theorem]{Definicja}
\newtheorem{lemma}[theorem]{Lemat}
\newtheorem{proposition}[theorem]{Stwierdzenie}
\newtheorem{corollary}[theorem]{Wniosek}
\newtheorem{remark}[theorem]{Uwaga}
\newtheorem{example}[theorem]{Przykład}
\fi

% --------------Dowód-------------------------------
\renewenvironment{proof}[1][\proofname]
{\par
  \vspace{-7pt}% remove the space after the theorem
  \pushQED{\qed}%
  \normalfont
  \topsep0pt \partopsep0pt % no space before
  \trivlist
  \item[\hskip\labelsep
        \sc % zastapic \itshape jak ktos woli napis "Dowód" kursywą
    #1\@addpunct{:}]\ignorespaces
}{%
  \popQED\endtrivlist\@endpefalse
  \addvspace{15pt} % some space after
}
\renewcommand{\qedsymbol}{$\square$}
%\renewcommand{\qedhere}{\hfill \qedsymbol}
% -----------------------------------------------------



% ABSTRACTS
% jak ktoś chce dwa streszczenia na jednej stronie, należy poniżej usunąć \clearpage



\newenvironment{streszczenie}{
\cleardoublepage
	\thispagestyle{empty}
  \vspace*{\fill}
  \begin{center}%
  \begin{localsize}{12}
    \textbf{Streszczenie} \\[10pt]
\noindent \@titleaux
  \end{localsize}
   \end{center}
}{\noindent \textbf{Słowa kluczowe:}~\@klucze  \vfill \clearpage}

\renewenvironment{abstract}{
	\cleardoublepage
	\thispagestyle{empty}
  \vspace*{\fill}
  \begin{center}%
  \begin{localsize}{12}
    \textbf{Abstract} \\[10pt]
\noindent \@title 
  \end{localsize}
  \end{center}}%
  {\noindent \textbf{Keywords:}~\@keywords  \vfill \clearpage}
 



%------------------------------------------------------------------------------

\ifmyclass@en
	\newcommand{\makestatement}
	{
	\cleardoublepage
	\thispagestyle{empty}
	%\pdfbookmark[0]{\@statementcap}{statement}

	\noindent\@author\hfill Warsaw, ..................\\
	Nr albumu \@album
	\par\vspace{5cm}
	\begin{center}
	  {Declaration}
	\end{center}

	%\par\vspace{10mm}
	I hereby declare that the thesis entitled ,,\@title '', submitted for the \@type\ degree, supervised  by \@supervisor , is entirely my original work apart from the recognized reference.
	\vspace{2cm}

	\begin{flushright}
 	 	\begin{minipage}{50mm}
   	 		\begin{center}
   	   		..............................................\\
    	 	 \vspace{5mm}
    	 	 \@author
   			 \end{center}
  		\end{minipage}
	\end{flushright}
	}
\else

	\newcommand{\makestatement}
	{
	\cleardoublepage
	\thispagestyle{empty}
%	\pdfbookmark[0]{\@statementcap}{statement}

	\noindent\@author\hfill Warszawa, dnia ..................\\
	Nr albumu \@album
	\par\vspace{5cm}
	\begin{center}
	  Oświadczenie
	\end{center}

	%\par\vspace{10mm}
	Oświadczam, że pracę \@type ką pod
	tytułem ,,\@titleaux '', której promotorem jest 	\@supervisor\
	wykonał{\ifthenelse{\boolean{lady}}{a}{e}}m
	samodzielnie, co poświadczam własnoręcznym podpisem.
	\vspace{2cm}

	\begin{flushright}
 	 	\begin{minipage}{50mm}
   	 		\begin{center}
   	   		..............................................\\
    	 	 \vspace{5mm}
    	 	 \@author
   			 \end{center}
  		\end{minipage}
	\end{flushright}
	}
\fi
%------------------------------------------------------------------------------%


\endinput
