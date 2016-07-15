#!/usr/bin/ruby
filename="./wiki.txt"
infile=File.new(filename,"r")
infile.each{
  |string1|
  string1.gsub!(/<br>/){""} #remove <br>
  string1.gsub!(/\\/){"\$\\backslash\$"}
  regex = Regexp.new(/'''.*?'''/) ########### Handeling bold ''' '''
  while (match = regex.match(string1)) != nil
    offsets = match.offset(0)
    startOfMatch = offsets[0]
    endOfMatch = offsets[1]
    string1[startOfMatch...endOfMatch] = "\\tb{" + match[0][3...match[0].length-3] + "}"
  end
  regex = Regexp.new(/''.*?''/) ########### Handeling italic '' ''
  while (match = regex.match(string1)) != nil
    offsets = match.offset(0)
    startOfMatch = offsets[0]
    endOfMatch = offsets[1]
    string1[startOfMatch...endOfMatch] = "\\textit{" + match[0][2...match[0].length-2] + "}"
  end
  regex = Regexp.new(/<tt>.*?<\/tt>/)
  while (match = regex.match(string1)) != nil
    offsets = match.offset(0)
    startOfMatch = offsets[0]
    endOfMatch = offsets[1]
    string1[startOfMatch...endOfMatch] = "\\ttt{" + match[0][4...match[0].length-5] + "}"
  end
  regex = Regexp.new(/===.*?===/) ########### Handeling subsections
  while (match = regex.match(string1)) != nil
    offsets = match.offset(0)
    startOfMatch = offsets[0]
    endOfMatch = offsets[1]
    string1[startOfMatch...endOfMatch] = "\\subsubsection{" + match[0][3...match[0].length-3] + "}"
  end
  regex = Regexp.new(/==.*?==/)
  while (match = regex.match(string1)) != nil
    offsets = match.offset(0)
    startOfMatch = offsets[0]
    endOfMatch = offsets[1]
    string1[startOfMatch...endOfMatch] = "\\subsection{" + match[0][2...match[0].length-2] + "}"
  end
  regex = Regexp.new(/=.*?=/)
  while (match = regex.match(string1)) != nil
    offsets = match.offset(0)
    startOfMatch = offsets[0]
    endOfMatch = offsets[1]
    string1[startOfMatch...endOfMatch] = "\\section{" + match[0][1...match[0].length-1] + "}"
  end
string1.gsub!(/([\#\$%_&~^{}])/){"\\" + $1} #### Handeling of special latex char # $ % & ~ _ { }, Warning ^ is incorrect
string1.gsub!(/\|/,'\textbar')
puts string1
}
