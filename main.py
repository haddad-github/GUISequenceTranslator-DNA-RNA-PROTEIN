import PySimpleGUI as sg

sg.theme("GreenMono")

layout = [

    #Enter sequence: [       ]
    [ sg.Text ("Enter DNA sequence:"), sg.Input(key = "seq_input")],

    #Reading direction [dropdown]
    [ sg.Text ("Reading direction"), sg.Combo( ["5'-3'", "3'-5'"], default_value= "5'-3'", key="inverted")],

    #ORF [dropdown]
    [ sg.Text ("ORF (which nuc. to start)"), sg.Combo( ["1", "2", "3"], default_value= "1", key="orf_input")],

    #(Translate)(Reset)
    [ sg.Button("Translate"), sg.Button("Reset") ],

    [ sg.Text ("RNA sequence:")],
    [ sg.Output(key="seq_output")],

    [ sg.Text ("Protein:")],
    [ sg.Output(key="prot_output")]


]

window = sg.Window('Sequence translator', layout)

while True:

    event, values = window.read()

    if event == sg.WINDOW_CLOSED:
        break

    #Reset all input and output values
    if event == "Reset":
        window['seq_input'].update("")
        window['seq_output'].update("")
        window['prot_output'].update("")

    if event == "Translate":

        #Protein equivalent for codons
        proteinDict = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
            'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
        }

        #RNA and DNA dictionaries
        rnaDict = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
        dnaDict = {'U': 'A', 'A': 'T', 'C': 'G', 'G': 'C'}

        #The only 4 accepted letters
        acceptedLetters = ['A', 'T', 'G', 'C']

        #Force it to upper case
        seqDNA = values['seq_input'].upper()

        #Translate DNA (not RNA) sequence into protein
        def getProtein(seq):

            protein = []

            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]

                if len(codon) == 3:
                    protein += proteinDict[codon]

            final_protein = "".join(protein)

            return final_protein

        #Clean unwanted letters
        for letter in seqDNA:
            if letter not in acceptedLetters:
                seqDNA = seqDNA.replace(letter, "")

        #Create RNA
        seqRNA = ""
        for letter in seqDNA:
            seqRNA += rnaDict[letter]

        #ORF (1st, 2nd or 3rd nucleotide as a starting point)
        #Nucleotide 1 = index 0; Nucleotide 2 = index 1; Nucleotide 3 = index 2
        combo2 = values['orf_input']

        if combo2 == '1':
            start = 0
            seqRNA = seqRNA[start::]

        if combo2 == '2':
            start = 1
            seqRNA = seqRNA[start::]

        if combo2 == '3':
            start = 2
            seqRNA = seqRNA[start::]

        #Inversion (3'-5')
        combo1 = values['inverted']

        if combo1 == "3'-5'":
            seqRNAreversed = ""
            for letter in seqRNA[::-1]:
                seqRNAreversed += letter
            window['seq_output'].update(seqRNAreversed)

            backToDNA = ""
            for letter in seqRNA:
                backToDNA += dnaDict[letter]

            backToDNA = backToDNA[::-1]

            myProt = getProtein(backToDNA)
            window['prot_output'].update(myProt)

        #Regular 5'-3'
        else:
            window['seq_output'].update(seqRNA)

            backToDNA = ""
            for letter in seqRNA:
                backToDNA += dnaDict[letter]

            myProt = getProtein(backToDNA)
            window['prot_output'].update(myProt)
