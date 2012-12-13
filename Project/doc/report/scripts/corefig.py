#!/usr/bin/env python

from __future__ import division

import colorsys

from collections import defaultdict


all_pos = {'L1': None,
'K1': None,
'J1': None,
'H1': None,
'G1': None,
'F1': None,
'E1': None,
'D1': None,
'C1': None,
'B1': None,
'A1': None,
'L2': None,
'K2': None,
'J2': None,
'H2': None,
'G2': None,
'F2': None,
'E2': None,
'D2': None,
'C2': None,
'B2': None,
'A2': None,
'L3': None,
'K3': None,
'J3': None,
'H3': None,
'G3': None,
'F3': None,
'E3': None,
'D3': None,
'C3': None,
'B3': None,
'A3': None,
'L4': None,
'K4': None,
'J4': None,
'H4': None,
'G4': None,
'F4': None,
'E4': None,
'D4': None,
'C4': None,
'B4': None,
'A4': None,
'L5': None,
'K5': None,
'J5': None,
'H5': None,
'G5': None,
'F5': None,
'E5': None,
'D5': None,
'C5': None,
'B5': None,
'A5': None,
'L6': None,
'K6': None,
'J6': None,
'H6': None,
'G6': None,
'F6': None,
'E6': None,
'D6': None,
'C6': None,
'B6': None,
'A6': None,
'L7': None,
'K7': None,
'J7': None,
'H7': None,
'G7': None,
'F7': None,
'E7': None,
'D7': None,
'C7': None,
'B7': None,
'A7': None,
'L8': None,
'K8': None,
'J8': None,
'H8': None,
'G8': None,
'F8': None,
'E8': None,
'D8': None,
'C8': None,
'B8': None,
'A8': None,
'L9': None,
'K9': None,
'J9': None,
'H9': None,
'G9': None,
'F9': None,
'E9': None,
'D9': None,
'C9': None,
'B9': None,
'A9': None,
'L10': None,
'K10': None,
'J10': None,
'H10': None,
'G10': None,
'F10': None,
'E10': None,
'D10': None,
'C10': None,
'B10': None,
'A10': None,
'L11': None,
'K11': None,
'J11': None,
'H11': None,
'G11': None,
'F11': None,
'E11': None,
'D11': None,
'C11': None,
'B11': None,
'A11': None,
}

numnames = {
    '1' : 'one',
    '2' : 'two',
    '3' : 'three',
    '4' : 'four',
    '5' : 'five',
    '6' : 'six',
    '7' : 'seven',
    '8' : 'eight',
    '9' : 'nine',
    '10' : 'ten',
    '11' : 'eleven',
    '12' : 'twelve',
    '13' : 'thirteen',
    '14' : 'fourteen',
    '15' : 'fifteen',
}

color_t = r"""\definecolor{{{colorname}}}{{rgb}}{{{r},{g},{b}}}
"""


main_fig_t = r"""
\begin{{figure}}[htbp]
    \centering
    
    % these dimensions are determined in arrow_dimms.ods

    \def\scale{{{scale}}}

    \def\latWidth{{0.2808363589*\scale}}
    
    \tikzset{{Assembly/.style={{
        inner sep=0pt,
        text width=\latWidth in,
        minimum size=\latWidth in,
        draw=black,
        align=center
        }}
    }}
    
    \def\matone{{blue}}
    \def\mattwo{{blue!50}}
    \def\matthree{{red}}
    \def\matfour{{yellow}}
    \def\matfive{{green}}
    \def\matsix{{yellow!50!red}}

    \scalebox{{{scalebox}}}{{

      \begin{{tikzpicture}}[x=1in,y=1in]
      
        % draw assembly row/column headers
        
        \draw[red, thick] ($(-5*\latWidth,7*\latWidth)$) node[above, anchor=south] {{A}} -- ($(-5*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-4*\latWidth,7*\latWidth)$) node[above, anchor=south] {{B}} -- ($(-4*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-3*\latWidth,7*\latWidth)$) node[above, anchor=south] {{C}} -- ($(-3*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-2*\latWidth,7*\latWidth)$) node[above, anchor=south] {{D}} -- ($(-2*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-1*\latWidth,7*\latWidth)$) node[above, anchor=south] {{E}} -- ($(-1*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(0*\latWidth,7*\latWidth)$) node[above, anchor=south] {{F}} -- ($(0*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(1*\latWidth,7*\latWidth)$) node[above, anchor=south] {{G}} -- ($(1*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(2*\latWidth,7*\latWidth)$) node[above, anchor=south] {{H}} -- ($(2*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(3*\latWidth,7*\latWidth)$) node[above, anchor=south] {{J}} -- ($(3*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(4*\latWidth,7*\latWidth)$) node[above, anchor=south] {{K}} -- ($(4*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(5*\latWidth,7*\latWidth)$) node[above, anchor=south] {{L}} -- ($(5*\latWidth,6*\latWidth)$);
        
        \begin{{scope}}[rotate=90]
          \draw[red, thick] ($(-5*\latWidth,7*\latWidth)$) node[left, anchor=east] {{1}} -- ($(-5*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-4*\latWidth,7*\latWidth)$) node[left, anchor=east] {{2}} -- ($(-4*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-3*\latWidth,7*\latWidth)$) node[left, anchor=east] {{3}} -- ($(-3*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-2*\latWidth,7*\latWidth)$) node[left, anchor=east] {{4}} -- ($(-2*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-1*\latWidth,7*\latWidth)$) node[left, anchor=east] {{5}} -- ($(-1*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(0*\latWidth,7*\latWidth)$) node[left, anchor=east] {{6}} -- ($(0*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(1*\latWidth,7*\latWidth)$) node[left, anchor=east] {{7}} -- ($(1*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(2*\latWidth,7*\latWidth)$) node[left, anchor=east] {{8}} -- ($(2*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(3*\latWidth,7*\latWidth)$) node[left, anchor=east] {{9}} -- ($(3*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(4*\latWidth,7*\latWidth)$) node[left, anchor=east] {{10}} -- ($(4*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(5*\latWidth,7*\latWidth)$) node[left, anchor=east] {{11}} -- ($(5*\latWidth,6*\latWidth)$);
        \end{{scope}}
        
        % draw fuel assembly nodes
        
        \node [Assembly, fill=\matfive]  at ($(-5*\latWidth,5*\latWidth)$)  {{{coredata[A11_content]}}}; % A11
        \node [Assembly, fill=\matfive]  at ($(-4*\latWidth,5*\latWidth)$)  {{{coredata[B11_content]}}}; % B11
        \node [Assembly, fill=\matfive]  at ($(-3*\latWidth,5*\latWidth)$)  {{{coredata[C11_content]}}}; % C11
        \node [Assembly, fill=\matfive]  at ($(-2*\latWidth,5*\latWidth)$)  {{{coredata[D11_content]}}}; % D11
        \node [Assembly, fill=\matfive]  at ($(-1*\latWidth,5*\latWidth)$)  {{{coredata[E11_content]}}}; % E11
        \node [Assembly, fill=\matfive]  at ($( 0*\latWidth,5*\latWidth)$)  {{{coredata[F11_content]}}}; % F11
        \node [Assembly, fill=\matfive]  at ($( 1*\latWidth,5*\latWidth)$)  {{{coredata[G11_content]}}}; % G11
        \node [Assembly, fill=\matfive]  at ($( 2*\latWidth,5*\latWidth)$)  {{{coredata[H11_content]}}}; % H11
        \node [Assembly, fill=\matfive]  at ($( 3*\latWidth,5*\latWidth)$)  {{{coredata[J11_content]}}}; % J11
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,5*\latWidth)$)  {{{coredata[K11_content]}}}; % K11
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,5*\latWidth)$)  {{{coredata[L11_content]}}}; % L11

        \node [Assembly, fill=\matfive]  at ($(-5*\latWidth,4*\latWidth)$)  {{{coredata[A10_content]}}}; % A10
        \node [Assembly, fill=\matfive]  at ($(-4*\latWidth,4*\latWidth)$)  {{{coredata[B10_content]}}}; % B10
        \node [Assembly, fill=\matfive]  at ($(-3*\latWidth,4*\latWidth)$)  {{{coredata[C10_content]}}}; % C10
        \node [Assembly, fill=\matfive]  at ($(-2*\latWidth,4*\latWidth)$)  {{{coredata[D10_content]}}}; % D10
        \node [Assembly, fill=\matfive]  at ($(-1*\latWidth,4*\latWidth)$)  {{{coredata[E10_content]}}}; % E10
        \node [Assembly, fill=\matfive]  at ($( 0*\latWidth,4*\latWidth)$)  {{{coredata[F10_content]}}}; % F10
        \node [Assembly, fill=\matfive]  at ($( 1*\latWidth,4*\latWidth)$)  {{{coredata[G10_content]}}}; % G10
        \node [Assembly, fill=\matfive]  at ($( 2*\latWidth,4*\latWidth)$)  {{{coredata[H10_content]}}}; % H10
        \node [Assembly, fill=\matfive]  at ($( 3*\latWidth,4*\latWidth)$)  {{{coredata[J10_content]}}}; % J10
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,4*\latWidth)$)  {{{coredata[K10_content]}}}; % K10
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,4*\latWidth)$)  {{{coredata[L10_content]}}}; % L10

        \node [Assembly, fill=\matthree] at ($(-5*\latWidth,3*\latWidth)$)  {{{coredata[A9_content]}}};  % A9
        \node [Assembly, fill=\matthree] at ($(-4*\latWidth,3*\latWidth)$)  {{{coredata[B9_content]}}};  % B9
        \node [Assembly, fill=\matthree] at ($(-3*\latWidth,3*\latWidth)$)  {{{coredata[C9_content]}}};  % C9
        \node [Assembly, fill=\matthree] at ($(-2*\latWidth,3*\latWidth)$)  {{{coredata[D9_content]}}};  % D9
        \node [Assembly, fill=\matthree] at ($(-1*\latWidth,3*\latWidth)$)  {{{coredata[E9_content]}}};  % E9
        \node [Assembly, fill=\matthree] at ($( 0*\latWidth,3*\latWidth)$)  {{{coredata[F9_content]}}};  % F9
        \node [Assembly, fill=\matthree] at ($( 1*\latWidth,3*\latWidth)$)  {{{coredata[G9_content]}}};  % G9
        \node [Assembly, fill=\matfive]  at ($( 2*\latWidth,3*\latWidth)$)  {{{coredata[H9_content]}}};  % H9
        \node [Assembly, fill=\matfive]  at ($( 3*\latWidth,3*\latWidth)$)  {{{coredata[J9_content]}}};  % J9
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,3*\latWidth)$)  {{{coredata[K9_content]}}};  % K9
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,3*\latWidth)$)  {{{coredata[L9_content]}}};  % L9

        \node [Assembly, fill=\matthree] at ($(-5*\latWidth,2*\latWidth)$)  {{{coredata[A8_content]}}};  % A8
        \node [Assembly, fill=\matthree] at ($(-4*\latWidth,2*\latWidth)$)  {{{coredata[B8_content]}}};  % B8
        \node [Assembly, fill=\matthree] at ($(-3*\latWidth,2*\latWidth)$)  {{{coredata[C8_content]}}};  % C8
        \node [Assembly, fill=\matthree] at ($(-2*\latWidth,2*\latWidth)$)  {{{coredata[D8_content]}}};  % D8
        \node [Assembly, fill=\matthree] at ($(-1*\latWidth,2*\latWidth)$)  {{{coredata[E8_content]}}};  % E8
        \node [Assembly, fill=\matthree] at ($( 0*\latWidth,2*\latWidth)$)  {{{coredata[F8_content]}}};  % F8
        \node [Assembly, fill=\matthree] at ($( 1*\latWidth,2*\latWidth)$)  {{{coredata[G8_content]}}};  % G8
        \node [Assembly, fill=\matfour]  at ($( 2*\latWidth,2*\latWidth)$)  {{{coredata[H8_content]}}};  % H8
        \node [Assembly, fill=\matfive]  at ($( 3*\latWidth,2*\latWidth)$)  {{{coredata[J8_content]}}};  % J8
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,2*\latWidth)$)  {{{coredata[K8_content]}}};  % K8
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,2*\latWidth)$)  {{{coredata[L8_content]}}};  % L8

        \node [Assembly, fill=\mattwo]   at ($(-5*\latWidth,1*\latWidth)$)  {{{coredata[A7_content]}}};  % A7
        \node [Assembly, fill=\matone]   at ($(-4*\latWidth,1*\latWidth)$)  {{{coredata[B7_content]}}};  % B7
        \node [Assembly, fill=\matone]   at ($(-3*\latWidth,1*\latWidth)$)  {{{coredata[C7_content]}}};  % C7
        \node [Assembly, fill=\matone]   at ($(-2*\latWidth,1*\latWidth)$)  {{{coredata[D7_content]}}};  % D7
        \node [Assembly, fill=\matone]   at ($(-1*\latWidth,1*\latWidth)$)  {{{coredata[E7_content]}}};  % E7
        \node [Assembly, fill=\mattwo]   at ($( 0*\latWidth,1*\latWidth)$)  {{{coredata[F7_content]}}};  % F7
        \node [Assembly, fill=\mattwo]   at ($( 1*\latWidth,1*\latWidth)$)  {{{coredata[G7_content]}}};  % G7
        \node [Assembly, fill=\matsix]   at ($( 2*\latWidth,1*\latWidth)$)  {{{coredata[H7_content]}}};  % H7
        \node [Assembly, fill=\matsix]   at ($( 3*\latWidth,1*\latWidth)$)  {{{coredata[K7_content]}}};  % J7
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,1*\latWidth)$)  {{{coredata[K7_content]}}};  % K7
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,1*\latWidth)$)  {{{coredata[L7_content]}}};  % L7

        \node [Assembly, fill=\mattwo]   at ($(-5*\latWidth,0*\latWidth)$)  {{{coredata[A6_content]}}};  % A6
        \node [Assembly, fill=\matone]   at ($(-4*\latWidth,0*\latWidth)$)  {{{coredata[B6_content]}}};  % B6
        \node [Assembly, fill=\matone]   at ($(-3*\latWidth,0*\latWidth)$)  {{{coredata[C6_content]}}};  % C6
        \node [Assembly, fill=\matone]   at ($(-2*\latWidth,0*\latWidth)$)  {{{coredata[D6_content]}}};  % D6
        \node [Assembly, fill=\matone]   at ($(-1*\latWidth,0*\latWidth)$)  {{{coredata[E6_content]}}};  % E6
        \node [Assembly, fill=\mattwo]   at ($( 0*\latWidth,0*\latWidth)$)  {{{coredata[F6_content]}}};  % F6
        \node [Assembly, fill=\mattwo]   at ($( 1*\latWidth,0*\latWidth)$)  {{{coredata[G6_content]}}};  % G6
        \node [Assembly, fill=\matsix]   at ($( 2*\latWidth,0*\latWidth)$)  {{{coredata[H6_content]}}};  % H6
        \node [Assembly, fill=\matsix]   at ($( 3*\latWidth,0*\latWidth)$)  {{{coredata[J6_content]}}};  % J6
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,0*\latWidth)$)  {{{coredata[K6_content]}}};  % K6
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,0*\latWidth)$)  {{{coredata[L6_content]}}};  % L6

        \node [Assembly, fill=\matone]   at ($(-5*\latWidth,-1*\latWidth)$) {{{coredata[A5_content]}}};  % A5
        \node [Assembly, fill=\matone]   at ($(-4*\latWidth,-1*\latWidth)$) {{{coredata[B5_content]}}};  % B5
        \node [Assembly, fill=\matone]   at ($(-3*\latWidth,-1*\latWidth)$) {{{coredata[C5_content]}}};  % C5
        \node [Assembly, fill=\matone]   at ($(-2*\latWidth,-1*\latWidth)$) {{{coredata[D5_content]}}};  % D5
        \node [Assembly, fill=\matone]   at ($(-1*\latWidth,-1*\latWidth)$) {{{coredata[E5_content]}}};  % E5
        \node [Assembly, fill=\matone]   at ($( 0*\latWidth,-1*\latWidth)$) {{{coredata[F5_content]}}};  % F5
        \node [Assembly, fill=\matone]   at ($( 1*\latWidth,-1*\latWidth)$) {{{coredata[G5_content]}}};  % G5
        \node [Assembly, fill=\matthree] at ($( 2*\latWidth,-1*\latWidth)$) {{{coredata[H5_content]}}};  % H5
        \node [Assembly, fill=\matthree] at ($( 3*\latWidth,-1*\latWidth)$) {{{coredata[J5_content]}}};  % J5
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,-1*\latWidth)$) {{{coredata[K5_content]}}};  % K5
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,-1*\latWidth)$) {{{coredata[L5_content]}}};  % L5

        \node [Assembly, fill=\matone]   at ($(-5*\latWidth,-2*\latWidth)$) {{{coredata[A4_content]}}};  % A4
        \node [Assembly, fill=\matone]   at ($(-4*\latWidth,-2*\latWidth)$) {{{coredata[B4_content]}}};  % B4
        \node [Assembly, fill=\matone]   at ($(-3*\latWidth,-2*\latWidth)$) {{{coredata[C4_content]}}};  % C4
        \node [Assembly, fill=\matone]   at ($(-2*\latWidth,-2*\latWidth)$) {{{coredata[D4_content]}}};  % D4
        \node [Assembly, fill=\matone]   at ($(-1*\latWidth,-2*\latWidth)$) {{{coredata[E4_content]}}};  % E4
        \node [Assembly, fill=\matone]   at ($( 0*\latWidth,-2*\latWidth)$) {{{coredata[F4_content]}}};  % F4
        \node [Assembly, fill=\matone]   at ($( 1*\latWidth,-2*\latWidth)$) {{{coredata[G4_content]}}};  % G4
        \node [Assembly, fill=\matthree] at ($( 2*\latWidth,-2*\latWidth)$) {{{coredata[H4_content]}}};  % H4
        \node [Assembly, fill=\matthree] at ($( 3*\latWidth,-2*\latWidth)$) {{{coredata[J4_content]}}};  % J4
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,-2*\latWidth)$) {{{coredata[K4_content]}}};  % K4
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,-2*\latWidth)$) {{{coredata[L4_content]}}};  % L4

        \node [Assembly, fill=\matone]   at ($(-5*\latWidth,-3*\latWidth)$) {{{coredata[A3_content]}}};  % A3
        \node [Assembly, fill=\matone]   at ($(-4*\latWidth,-3*\latWidth)$) {{{coredata[B3_content]}}};  % B3
        \node [Assembly, fill=\matone]   at ($(-3*\latWidth,-3*\latWidth)$) {{{coredata[C3_content]}}};  % C3
        \node [Assembly, fill=\matone]   at ($(-2*\latWidth,-3*\latWidth)$) {{{coredata[D3_content]}}};  % D3
        \node [Assembly, fill=\matone]   at ($(-1*\latWidth,-3*\latWidth)$) {{{coredata[E3_content]}}};  % E3
        \node [Assembly, fill=\matone]   at ($( 0*\latWidth,-3*\latWidth)$) {{{coredata[F3_content]}}};  % F3
        \node [Assembly, fill=\matone]   at ($( 1*\latWidth,-3*\latWidth)$) {{{coredata[G3_content]}}};  % G3
        \node [Assembly, fill=\matthree] at ($( 2*\latWidth,-3*\latWidth)$) {{{coredata[H3_content]}}};  % H3
        \node [Assembly, fill=\matthree] at ($( 3*\latWidth,-3*\latWidth)$) {{{coredata[J3_content]}}};  % J3
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,-3*\latWidth)$) {{{coredata[K3_content]}}};  % K3
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,-3*\latWidth)$) {{{coredata[L3_content]}}};  % L3

        \node [Assembly, fill=\matone]   at ($(-5*\latWidth,-4*\latWidth)$) {{{coredata[A2_content]}}};  % A2
        \node [Assembly, fill=\matone]   at ($(-4*\latWidth,-4*\latWidth)$) {{{coredata[B2_content]}}};  % B2
        \node [Assembly, fill=\matone]   at ($(-3*\latWidth,-4*\latWidth)$) {{{coredata[C2_content]}}};  % C2
        \node [Assembly, fill=\matone]   at ($(-2*\latWidth,-4*\latWidth)$) {{{coredata[D2_content]}}};  % D2
        \node [Assembly, fill=\matone]   at ($(-1*\latWidth,-4*\latWidth)$) {{{coredata[E2_content]}}};  % E2
        \node [Assembly, fill=\matone]   at ($( 0*\latWidth,-4*\latWidth)$) {{{coredata[F2_content]}}};  % F2
        \node [Assembly, fill=\matone]   at ($( 1*\latWidth,-4*\latWidth)$) {{{coredata[G2_content]}}};  % G2
        \node [Assembly, fill=\matthree] at ($( 2*\latWidth,-4*\latWidth)$) {{{coredata[H2_content]}}};  % H2
        \node [Assembly, fill=\matthree] at ($( 3*\latWidth,-4*\latWidth)$) {{{coredata[J2_content]}}};  % J2
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,-4*\latWidth)$) {{{coredata[K2_content]}}};  % K2
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,-4*\latWidth)$) {{{coredata[L2_content]}}};  % L2

        \node [Assembly, fill=\mattwo]   at ($(-5*\latWidth,-5*\latWidth)$) {{{coredata[A1_content]}}};  % A1
        \node [Assembly, fill=\matone]   at ($(-4*\latWidth,-5*\latWidth)$) {{{coredata[B1_content]}}};  % B1
        \node [Assembly, fill=\matone]   at ($(-3*\latWidth,-5*\latWidth)$) {{{coredata[C1_content]}}};  % C1
        \node [Assembly, fill=\matone]   at ($(-2*\latWidth,-5*\latWidth)$) {{{coredata[D1_content]}}};  % D1
        \node [Assembly, fill=\matone]   at ($(-1*\latWidth,-5*\latWidth)$) {{{coredata[E1_content]}}};  % E1
        \node [Assembly, fill=\mattwo]   at ($( 0*\latWidth,-5*\latWidth)$) {{{coredata[F1_content]}}};  % F1
        \node [Assembly, fill=\mattwo]   at ($( 1*\latWidth,-5*\latWidth)$) {{{coredata[G1_content]}}};  % G1
        \node [Assembly, fill=\matthree] at ($( 2*\latWidth,-5*\latWidth)$) {{{coredata[H1_content]}}};  % H1
        \node [Assembly, fill=\matthree] at ($( 3*\latWidth,-5*\latWidth)$) {{{coredata[J1_content]}}};  % J1
        \node [Assembly, fill=\matfive]  at ($( 4*\latWidth,-5*\latWidth)$) {{{coredata[K1_content]}}};  % K1
        \node [Assembly, fill=\matfive]  at ($( 5*\latWidth,-5*\latWidth)$) {{{coredata[L1_content]}}}; `% L1

      \end{{tikzpicture}}
    }}
    
{legend}

    \caption{altcap}{{{caption} \label{{{label}}}}}
\end{{figure}}

"""

main_fig_clr_t = r"""
\begin{{figure}}[htbp]
    \centering
    
    % these dimensions are determined in arrow_dimms.ods

    \def\scale{{{scale}}}

    \def\latWidth{{0.2808363589*\scale}}
    
    \tikzset{{Assembly/.style={{
        inner sep=0pt,
        text width=\latWidth in,
        minimum size=\latWidth in,
        draw=black,
        align=center
        }}
    }}
    
    \def\matone{{blue}}
    \def\mattwo{{blue!50}}
    \def\matthree{{red}}
    \def\matfour{{yellow}}
    \def\matfive{{green}}
    \def\matsix{{yellow!50!red}}

    \scalebox{{{scalebox}}}{{

      \begin{{tikzpicture}}[x=1in,y=1in]
      
        % draw assembly row/column headers
        
        \draw[red, thick] ($(-5*\latWidth,7*\latWidth)$) node[above, anchor=south] {{A}} -- ($(-5*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-4*\latWidth,7*\latWidth)$) node[above, anchor=south] {{B}} -- ($(-4*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-3*\latWidth,7*\latWidth)$) node[above, anchor=south] {{C}} -- ($(-3*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-2*\latWidth,7*\latWidth)$) node[above, anchor=south] {{D}} -- ($(-2*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(-1*\latWidth,7*\latWidth)$) node[above, anchor=south] {{E}} -- ($(-1*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(0*\latWidth,7*\latWidth)$) node[above, anchor=south] {{F}} -- ($(0*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(1*\latWidth,7*\latWidth)$) node[above, anchor=south] {{G}} -- ($(1*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(2*\latWidth,7*\latWidth)$) node[above, anchor=south] {{H}} -- ($(2*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(3*\latWidth,7*\latWidth)$) node[above, anchor=south] {{J}} -- ($(3*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(4*\latWidth,7*\latWidth)$) node[above, anchor=south] {{K}} -- ($(4*\latWidth,6*\latWidth)$);
        \draw[red, thick] ($(5*\latWidth,7*\latWidth)$) node[above, anchor=south] {{L}} -- ($(5*\latWidth,6*\latWidth)$);
        
        \begin{{scope}}[rotate=90]
          \draw[red, thick] ($(-5*\latWidth,7*\latWidth)$) node[left, anchor=east] {{1}} -- ($(-5*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-4*\latWidth,7*\latWidth)$) node[left, anchor=east] {{2}} -- ($(-4*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-3*\latWidth,7*\latWidth)$) node[left, anchor=east] {{3}} -- ($(-3*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-2*\latWidth,7*\latWidth)$) node[left, anchor=east] {{4}} -- ($(-2*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(-1*\latWidth,7*\latWidth)$) node[left, anchor=east] {{5}} -- ($(-1*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(0*\latWidth,7*\latWidth)$) node[left, anchor=east] {{6}} -- ($(0*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(1*\latWidth,7*\latWidth)$) node[left, anchor=east] {{7}} -- ($(1*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(2*\latWidth,7*\latWidth)$) node[left, anchor=east] {{8}} -- ($(2*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(3*\latWidth,7*\latWidth)$) node[left, anchor=east] {{9}} -- ($(3*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(4*\latWidth,7*\latWidth)$) node[left, anchor=east] {{10}} -- ($(4*\latWidth,6*\latWidth)$);
          \draw[red, thick] ($(5*\latWidth,7*\latWidth)$) node[left, anchor=east] {{11}} -- ($(5*\latWidth,6*\latWidth)$);
        \end{{scope}}
        
        % draw fuel assembly nodes
        
        \node [Assembly, fill=Aelevencolor]  at ($(-5*\latWidth,5*\latWidth)$)  {{{coredata[A11_content]}}}; % A11
        \node [Assembly, fill=Belevencolor]  at ($(-4*\latWidth,5*\latWidth)$)  {{{coredata[B11_content]}}}; % B11
        \node [Assembly, fill=Celevencolor]  at ($(-3*\latWidth,5*\latWidth)$)  {{{coredata[C11_content]}}}; % C11
        \node [Assembly, fill=Delevencolor]  at ($(-2*\latWidth,5*\latWidth)$)  {{{coredata[D11_content]}}}; % D11
        \node [Assembly, fill=Eelevencolor]  at ($(-1*\latWidth,5*\latWidth)$)  {{{coredata[E11_content]}}}; % E11
        \node [Assembly, fill=Felevencolor]  at ($( 0*\latWidth,5*\latWidth)$)  {{{coredata[F11_content]}}}; % F11
        \node [Assembly, fill=Gelevencolor]  at ($( 1*\latWidth,5*\latWidth)$)  {{{coredata[G11_content]}}}; % G11
        \node [Assembly, fill=Helevencolor]  at ($( 2*\latWidth,5*\latWidth)$)  {{{coredata[H11_content]}}}; % H11
        \node [Assembly, fill=Jelevencolor]  at ($( 3*\latWidth,5*\latWidth)$)  {{{coredata[J11_content]}}}; % J11
        \node [Assembly, fill=Kelevencolor]  at ($( 4*\latWidth,5*\latWidth)$)  {{{coredata[K11_content]}}}; % K11
        \node [Assembly, fill=Lelevencolor]  at ($( 5*\latWidth,5*\latWidth)$)  {{{coredata[L11_content]}}}; % L11

        \node [Assembly, fill=Atencolor]  at ($(-5*\latWidth,4*\latWidth)$)  {{{coredata[A10_content]}}}; % A10
        \node [Assembly, fill=Btencolor]  at ($(-4*\latWidth,4*\latWidth)$)  {{{coredata[B10_content]}}}; % B10
        \node [Assembly, fill=Ctencolor]  at ($(-3*\latWidth,4*\latWidth)$)  {{{coredata[C10_content]}}}; % C10
        \node [Assembly, fill=Dtencolor]  at ($(-2*\latWidth,4*\latWidth)$)  {{{coredata[D10_content]}}}; % D10
        \node [Assembly, fill=Etencolor]  at ($(-1*\latWidth,4*\latWidth)$)  {{{coredata[E10_content]}}}; % E10
        \node [Assembly, fill=Ftencolor]  at ($( 0*\latWidth,4*\latWidth)$)  {{{coredata[F10_content]}}}; % F10
        \node [Assembly, fill=Gtencolor]  at ($( 1*\latWidth,4*\latWidth)$)  {{{coredata[G10_content]}}}; % G10
        \node [Assembly, fill=Htencolor]  at ($( 2*\latWidth,4*\latWidth)$)  {{{coredata[H10_content]}}}; % H10
        \node [Assembly, fill=Jtencolor]  at ($( 3*\latWidth,4*\latWidth)$)  {{{coredata[J10_content]}}}; % J10
        \node [Assembly, fill=Ktencolor]  at ($( 4*\latWidth,4*\latWidth)$)  {{{coredata[K10_content]}}}; % K10
        \node [Assembly, fill=Ltencolor]  at ($( 5*\latWidth,4*\latWidth)$)  {{{coredata[L10_content]}}}; % L10

        \node [Assembly, fill=Aninecolor] at ($(-5*\latWidth,3*\latWidth)$)  {{{coredata[A9_content]}}};  % A9
        \node [Assembly, fill=Bninecolor] at ($(-4*\latWidth,3*\latWidth)$)  {{{coredata[B9_content]}}};  % B9
        \node [Assembly, fill=Cninecolor] at ($(-3*\latWidth,3*\latWidth)$)  {{{coredata[C9_content]}}};  % C9
        \node [Assembly, fill=Dninecolor] at ($(-2*\latWidth,3*\latWidth)$)  {{{coredata[D9_content]}}};  % D9
        \node [Assembly, fill=Eninecolor] at ($(-1*\latWidth,3*\latWidth)$)  {{{coredata[E9_content]}}};  % E9
        \node [Assembly, fill=Fninecolor] at ($( 0*\latWidth,3*\latWidth)$)  {{{coredata[F9_content]}}};  % F9
        \node [Assembly, fill=Gninecolor] at ($( 1*\latWidth,3*\latWidth)$)  {{{coredata[G9_content]}}};  % G9
        \node [Assembly, fill=Hninecolor] at ($( 2*\latWidth,3*\latWidth)$)  {{{coredata[H9_content]}}};  % H9
        \node [Assembly, fill=Jninecolor] at ($( 3*\latWidth,3*\latWidth)$)  {{{coredata[J9_content]}}};  % J9
        \node [Assembly, fill=Kninecolor] at ($( 4*\latWidth,3*\latWidth)$)  {{{coredata[K9_content]}}};  % K9
        \node [Assembly, fill=Lninecolor] at ($( 5*\latWidth,3*\latWidth)$)  {{{coredata[L9_content]}}};  % L9

        \node [Assembly, fill=Aeightcolor] at ($(-5*\latWidth,2*\latWidth)$)  {{{coredata[A8_content]}}};  % A8
        \node [Assembly, fill=Beightcolor] at ($(-4*\latWidth,2*\latWidth)$)  {{{coredata[B8_content]}}};  % B8
        \node [Assembly, fill=Ceightcolor] at ($(-3*\latWidth,2*\latWidth)$)  {{{coredata[C8_content]}}};  % C8
        \node [Assembly, fill=Deightcolor] at ($(-2*\latWidth,2*\latWidth)$)  {{{coredata[D8_content]}}};  % D8
        \node [Assembly, fill=Eeightcolor] at ($(-1*\latWidth,2*\latWidth)$)  {{{coredata[E8_content]}}};  % E8
        \node [Assembly, fill=Feightcolor] at ($( 0*\latWidth,2*\latWidth)$)  {{{coredata[F8_content]}}};  % F8
        \node [Assembly, fill=Geightcolor] at ($( 1*\latWidth,2*\latWidth)$)  {{{coredata[G8_content]}}};  % G8
        \node [Assembly, fill=Heightcolor] at ($( 2*\latWidth,2*\latWidth)$)  {{{coredata[H8_content]}}};  % H8
        \node [Assembly, fill=Jeightcolor] at ($( 3*\latWidth,2*\latWidth)$)  {{{coredata[J8_content]}}};  % J8
        \node [Assembly, fill=Keightcolor] at ($( 4*\latWidth,2*\latWidth)$)  {{{coredata[K8_content]}}};  % K8
        \node [Assembly, fill=Leightcolor] at ($( 5*\latWidth,2*\latWidth)$)  {{{coredata[L8_content]}}};  % L8

        \node [Assembly, fill=Asevencolor]  at ($(-5*\latWidth,1*\latWidth)$)  {{{coredata[A7_content]}}};  % A7
        \node [Assembly, fill=Bsevencolor]  at ($(-4*\latWidth,1*\latWidth)$)  {{{coredata[B7_content]}}};  % B7
        \node [Assembly, fill=Csevencolor]  at ($(-3*\latWidth,1*\latWidth)$)  {{{coredata[C7_content]}}};  % C7
        \node [Assembly, fill=Dsevencolor]  at ($(-2*\latWidth,1*\latWidth)$)  {{{coredata[D7_content]}}};  % D7
        \node [Assembly, fill=Esevencolor]  at ($(-1*\latWidth,1*\latWidth)$)  {{{coredata[E7_content]}}};  % E7
        \node [Assembly, fill=Fsevencolor]  at ($( 0*\latWidth,1*\latWidth)$)  {{{coredata[F7_content]}}};  % F7
        \node [Assembly, fill=Gsevencolor]  at ($( 1*\latWidth,1*\latWidth)$)  {{{coredata[G7_content]}}};  % G7
        \node [Assembly, fill=Hsevencolor]  at ($( 2*\latWidth,1*\latWidth)$)  {{{coredata[H7_content]}}};  % H7
        \node [Assembly, fill=Jsevencolor]  at ($( 3*\latWidth,1*\latWidth)$)  {{{coredata[J7_content]}}};  % J7
        \node [Assembly, fill=Ksevencolor]  at ($( 4*\latWidth,1*\latWidth)$)  {{{coredata[K7_content]}}};  % K7
        \node [Assembly, fill=Lsevencolor]  at ($( 5*\latWidth,1*\latWidth)$)  {{{coredata[L7_content]}}};  % L7

        \node [Assembly, fill=Asixcolor]  at ($(-5*\latWidth,0*\latWidth)$)  {{{coredata[A6_content]}}};  % A6
        \node [Assembly, fill=Bsixcolor]  at ($(-4*\latWidth,0*\latWidth)$)  {{{coredata[B6_content]}}};  % B6
        \node [Assembly, fill=Csixcolor]  at ($(-3*\latWidth,0*\latWidth)$)  {{{coredata[C6_content]}}};  % C6
        \node [Assembly, fill=Dsixcolor]  at ($(-2*\latWidth,0*\latWidth)$)  {{{coredata[D6_content]}}};  % D6
        \node [Assembly, fill=Esixcolor]  at ($(-1*\latWidth,0*\latWidth)$)  {{{coredata[E6_content]}}};  % E6
        \node [Assembly, fill=Fsixcolor]  at ($( 0*\latWidth,0*\latWidth)$)  {{{coredata[F6_content]}}};  % F6
        \node [Assembly, fill=Gsixcolor]  at ($( 1*\latWidth,0*\latWidth)$)  {{{coredata[G6_content]}}};  % G6
        \node [Assembly, fill=Hsixcolor]  at ($( 2*\latWidth,0*\latWidth)$)  {{{coredata[H6_content]}}};  % H6
        \node [Assembly, fill=Jsixcolor]  at ($( 3*\latWidth,0*\latWidth)$)  {{{coredata[J6_content]}}};  % J6
        \node [Assembly, fill=Ksixcolor]  at ($( 4*\latWidth,0*\latWidth)$)  {{{coredata[K6_content]}}};  % K6
        \node [Assembly, fill=Lsixcolor]  at ($( 5*\latWidth,0*\latWidth)$)  {{{coredata[L6_content]}}};  % L6

        \node [Assembly, fill=Afivecolor]   at ($(-5*\latWidth,-1*\latWidth)$) {{{coredata[A5_content]}}};  % A5
        \node [Assembly, fill=Bfivecolor]   at ($(-4*\latWidth,-1*\latWidth)$) {{{coredata[B5_content]}}};  % B5
        \node [Assembly, fill=Cfivecolor]   at ($(-3*\latWidth,-1*\latWidth)$) {{{coredata[C5_content]}}};  % C5
        \node [Assembly, fill=Dfivecolor]   at ($(-2*\latWidth,-1*\latWidth)$) {{{coredata[D5_content]}}};  % D5
        \node [Assembly, fill=Efivecolor]   at ($(-1*\latWidth,-1*\latWidth)$) {{{coredata[E5_content]}}};  % E5
        \node [Assembly, fill=Ffivecolor]   at ($( 0*\latWidth,-1*\latWidth)$) {{{coredata[F5_content]}}};  % F5
        \node [Assembly, fill=Gfivecolor]   at ($( 1*\latWidth,-1*\latWidth)$) {{{coredata[G5_content]}}};  % G5
        \node [Assembly, fill=Hfivecolor] at ($( 2*\latWidth,-1*\latWidth)$) {{{coredata[H5_content]}}};  % H5
        \node [Assembly, fill=Jfivecolor] at ($( 3*\latWidth,-1*\latWidth)$) {{{coredata[J5_content]}}};  % J5
        \node [Assembly, fill=Kfivecolor]  at ($( 4*\latWidth,-1*\latWidth)$) {{{coredata[K5_content]}}};  % K5
        \node [Assembly, fill=Lfivecolor]  at ($( 5*\latWidth,-1*\latWidth)$) {{{coredata[L5_content]}}};  % L5

        \node [Assembly, fill=Afourcolor]   at ($(-5*\latWidth,-2*\latWidth)$) {{{coredata[A4_content]}}};  % A4
        \node [Assembly, fill=Bfourcolor]   at ($(-4*\latWidth,-2*\latWidth)$) {{{coredata[B4_content]}}};  % B4
        \node [Assembly, fill=Cfourcolor]   at ($(-3*\latWidth,-2*\latWidth)$) {{{coredata[C4_content]}}};  % C4
        \node [Assembly, fill=Dfourcolor]   at ($(-2*\latWidth,-2*\latWidth)$) {{{coredata[D4_content]}}};  % D4
        \node [Assembly, fill=Efourcolor]   at ($(-1*\latWidth,-2*\latWidth)$) {{{coredata[E4_content]}}};  % E4
        \node [Assembly, fill=Ffourcolor]   at ($( 0*\latWidth,-2*\latWidth)$) {{{coredata[F4_content]}}};  % F4
        \node [Assembly, fill=Gfourcolor]   at ($( 1*\latWidth,-2*\latWidth)$) {{{coredata[G4_content]}}};  % G4
        \node [Assembly, fill=Hfourcolor] at ($( 2*\latWidth,-2*\latWidth)$) {{{coredata[H4_content]}}};  % H4
        \node [Assembly, fill=Jfourcolor] at ($( 3*\latWidth,-2*\latWidth)$) {{{coredata[J4_content]}}};  % J4
        \node [Assembly, fill=Kfourcolor]  at ($( 4*\latWidth,-2*\latWidth)$) {{{coredata[K4_content]}}};  % K4
        \node [Assembly, fill=Lfourcolor]  at ($( 5*\latWidth,-2*\latWidth)$) {{{coredata[L4_content]}}};  % L4

        \node [Assembly, fill=Athreecolor]   at ($(-5*\latWidth,-3*\latWidth)$) {{{coredata[A3_content]}}};  % A3
        \node [Assembly, fill=Bthreecolor]   at ($(-4*\latWidth,-3*\latWidth)$) {{{coredata[B3_content]}}};  % B3
        \node [Assembly, fill=Cthreecolor]   at ($(-3*\latWidth,-3*\latWidth)$) {{{coredata[C3_content]}}};  % C3
        \node [Assembly, fill=Dthreecolor]   at ($(-2*\latWidth,-3*\latWidth)$) {{{coredata[D3_content]}}};  % D3
        \node [Assembly, fill=Ethreecolor]   at ($(-1*\latWidth,-3*\latWidth)$) {{{coredata[E3_content]}}};  % E3
        \node [Assembly, fill=Fthreecolor]   at ($( 0*\latWidth,-3*\latWidth)$) {{{coredata[F3_content]}}};  % F3
        \node [Assembly, fill=Gthreecolor]   at ($( 1*\latWidth,-3*\latWidth)$) {{{coredata[G3_content]}}};  % G3
        \node [Assembly, fill=Hthreecolor] at ($( 2*\latWidth,-3*\latWidth)$) {{{coredata[H3_content]}}};  % H3
        \node [Assembly, fill=Jthreecolor] at ($( 3*\latWidth,-3*\latWidth)$) {{{coredata[J3_content]}}};  % J3
        \node [Assembly, fill=Kthreecolor]  at ($( 4*\latWidth,-3*\latWidth)$) {{{coredata[K3_content]}}};  % K3
        \node [Assembly, fill=Lthreecolor]  at ($( 5*\latWidth,-3*\latWidth)$) {{{coredata[L3_content]}}};  % L3

        \node [Assembly, fill=Atwocolor]   at ($(-5*\latWidth,-4*\latWidth)$) {{{coredata[A2_content]}}};  % A2
        \node [Assembly, fill=Btwocolor]   at ($(-4*\latWidth,-4*\latWidth)$) {{{coredata[B2_content]}}};  % B2
        \node [Assembly, fill=Ctwocolor]   at ($(-3*\latWidth,-4*\latWidth)$) {{{coredata[C2_content]}}};  % C2
        \node [Assembly, fill=Dtwocolor]   at ($(-2*\latWidth,-4*\latWidth)$) {{{coredata[D2_content]}}};  % D2
        \node [Assembly, fill=Etwocolor]   at ($(-1*\latWidth,-4*\latWidth)$) {{{coredata[E2_content]}}};  % E2
        \node [Assembly, fill=Ftwocolor]   at ($( 0*\latWidth,-4*\latWidth)$) {{{coredata[F2_content]}}};  % F2
        \node [Assembly, fill=Gtwocolor]   at ($( 1*\latWidth,-4*\latWidth)$) {{{coredata[G2_content]}}};  % G2
        \node [Assembly, fill=Htwocolor] at ($( 2*\latWidth,-4*\latWidth)$) {{{coredata[H2_content]}}};  % H2
        \node [Assembly, fill=Jtwocolor] at ($( 3*\latWidth,-4*\latWidth)$) {{{coredata[J2_content]}}};  % J2
        \node [Assembly, fill=Ktwocolor]  at ($( 4*\latWidth,-4*\latWidth)$) {{{coredata[K2_content]}}};  % K2
        \node [Assembly, fill=Ltwocolor]  at ($( 5*\latWidth,-4*\latWidth)$) {{{coredata[L2_content]}}};  % L2

        \node [Assembly, fill=Aonecolor]   at ($(-5*\latWidth,-5*\latWidth)$) {{{coredata[A1_content]}}};  % A1
        \node [Assembly, fill=Bonecolor]   at ($(-4*\latWidth,-5*\latWidth)$) {{{coredata[B1_content]}}};  % B1
        \node [Assembly, fill=Conecolor]   at ($(-3*\latWidth,-5*\latWidth)$) {{{coredata[C1_content]}}};  % C1
        \node [Assembly, fill=Donecolor]   at ($(-2*\latWidth,-5*\latWidth)$) {{{coredata[D1_content]}}};  % D1
        \node [Assembly, fill=Eonecolor]   at ($(-1*\latWidth,-5*\latWidth)$) {{{coredata[E1_content]}}};  % E1
        \node [Assembly, fill=Fonecolor]   at ($( 0*\latWidth,-5*\latWidth)$) {{{coredata[F1_content]}}};  % F1
        \node [Assembly, fill=Gonecolor]   at ($( 1*\latWidth,-5*\latWidth)$) {{{coredata[G1_content]}}};  % G1
        \node [Assembly, fill=Honecolor] at ($( 2*\latWidth,-5*\latWidth)$) {{{coredata[H1_content]}}};  % H1
        \node [Assembly, fill=Jonecolor] at ($( 3*\latWidth,-5*\latWidth)$) {{{coredata[J1_content]}}};  % J1
        \node [Assembly, fill=Konecolor]  at ($( 4*\latWidth,-5*\latWidth)$) {{{coredata[K1_content]}}};  % K1
        \node [Assembly, fill=Lonecolor]  at ($( 5*\latWidth,-5*\latWidth)$) {{{coredata[L1_content]}}}; `% L1
        
      \end{{tikzpicture}}
    }}
    
{legend}

    \caption{altcap}{{{caption} \label{{{label}}}}}
\end{{figure}}

"""



legend_t = r"""    % make the legend
    \begin{{tikzpicture}}
      \matrix [matrix of nodes]
          {{
              \node [Assembly, fill=\matone] at (0,0) {{}}; & Fuel 1, blade in~~~ & 
              \node [Assembly, fill=\mattwo] at (0,0) {{}}; & Fuel 1, blade out~~ &
              \node [Assembly, fill=\matthree] at (0,0) {{}}; & Fuel 2, blade in~~~ \\
              \node [Assembly, fill=\matfour] at (0,0) {{}}; & Fuel 2, blade out~~~ & 
              \node [Assembly, fill=\matsix] at (0,0) {{}}; & Fuel 2, blade in/out~~~&
              \node [Assembly, fill=\matfive] at (0,0) {{}}; & Reflector~~\\ 
          }};
    \end{{tikzpicture}}"""

hyper_t = r""", hyperlink node={target}"""


class CoreFig:
  def __init__(self,caption,label,altcap='',scale=1.0,scalebox=1.0,colorbyval=False):
    self.datadict = defaultdict(str)
    self.legend = ""
    self.caption = caption
    if altcap != '':
      self.altcap = '[{0}]'.format(altcap)
    else:
      self.altcap = ''
    self.label = label
    self.quarter = False
    self.structs = True
    self.scale=scale
    self.scalebox=scalebox
    self.colorbyval=colorbyval
    self.positions = all_pos
    self.values = []

  def set_quarter(self,sett):
    self.quarter = sett

  def set_structs(self,sett):
    self.structs = sett

  def set_legend(self):
    self.legend = legend_t.format()
    
  def clear_legend(self):
    self.legend = ""

  def set_pos(self,pos,tex,hyper=None,val=None):
    """Set the tex in the tikz node at assembly position pos.

      args
        pos - Core position, e.g.: 'E15'
        tex - The tex to put in the node, e.g.  1.2//3.4//5.6   or  $\mathrm{S}_\mathrm{A}$, etc

     kwargs
        hyper - Default to none - name of a hypertarget for this node to link to
    
    """
    if not pos in self.positions: raise Exception("Unknown position {0}".format(pos))
    
    self.datadict[pos+"_content"] = tex
    if hyper:
      self.datadict[pos+"_hyper"] = hyper_t.format(target=hyper)

    if self.colorbyval:
      self.positions[pos] = val
      if val:
        self.values.append(val)


  def write_fig(self,filename):
    if self.quarter:
      if self.colorbyval:
        templ = quart_fig_clr_t
      else:
        templ = quart_fig_t
    else:
      if self.colorbyval:
        templ = main_fig_clr_t
      else:
        templ = main_fig_t

    

    with open(filename,'w') as fh:

      if self.colorbyval:
        for pos,val in self.positions.items():
          p = pos.replace(pos[1:],numnames[pos[1:]])
          pocolkey = "{0}color".format(p)
          if not val:
            fh.write(color_t.format(colorname=pocolkey,r=1,g=1,b=1))
          else:
            h = (1-(val-min(self.values))/(max(self.values) - min(self.values)))/3
            s = 0.8
            v = 0.8
            r,g,b = colorsys.hsv_to_rgb(h,s,v)
            fh.write(color_t.format(colorname=pocolkey,r=r,g=g,b=b))
    
      head_lhm = r"\RPVOR/\latWidth"
      head_llm = "1"
      struct = ""
      
      if not self.structs:
        head_lhm = "1"
        head_llm = "0.6"
        struct = "%"
        
        if not self.quarter:
          head_lhm = "8.5"
    
      fh.write(templ.format(coredata=self.datadict,
                            legend=self.legend,
                            caption=self.caption,
                            altcap=self.altcap,
                            label=self.label,
                            scale=self.scale,
                            scalebox=self.scalebox,
                            head_lhm=head_lhm,
                            head_llm=head_llm,
                            struct=struct))

def main():

  fig = CoreFig('test caption','test label',scale=2.0,scalebox=0.8)
  fig.set_legend()
  fig.set_quarter(True)
  fig.set_pos('D11','test')
  fig.set_pos('D12','test2',hyper='mytarget')
  fig.write_fig('tmp')

if __name__ == "__main__":
  main()
