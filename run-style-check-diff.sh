find . -name "*.cpp" -o -name "*.hpp" | while read file; 
do clang-format -style=file -output-replacements-xml "$file" > "$file.xml" 
if grep -q '<replacement ' "$file.xml"; then echo "File $file not formatted correctly" 
exit 1 
fi 
done
