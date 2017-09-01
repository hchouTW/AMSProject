#!/usr/bin/ruby
def select(name, options)
	if (options.length == 0) 
		return ""
	end
	@id = 0
	puts "#{name} : "
	options.each.with_index { |x, i| puts "#{i+1}) #{x}"}
	while @id < 1 || @id > options.length
		print "#? "
		@id = $stdin.gets.chomp.to_i
	end
	puts "Select (#{@id}) #{options[@id-1]}"
	puts ""
	return options[@id-1]
end

def inputValue(mode, name, inval, region)
	print "#{name} : "
	if (mode.eql? "UI")
		@val = $stdin.gets.chomp.to_i
		while @val < region[0] || @val > region[1]
			print "#{name} : "
			@val = $stdin.gets.chomp.to_i
		end
		print @val
		puts ""
		return @val
	elsif (mode.eql? "CONFIG")
		if !(inval.to_i.to_s == inval.to_s)
			puts "{ #{inval} } isn't a integer."
			exit
		end
		@val = inval.to_i
		if (@val < region[0]) || (@val > region[1])
			puts "{ #{@val} } is out the region (#{region[0]} ... #{region[1]})"
			exit
		end
		print @val
		puts ""
		return @val
	end
	puts ""
	return 0
end

def show(name, option, options = "")
	if (options.length == 0)
		puts "#{name} : #{option}\n\n"
	else
		@isfind = 0
		options.each do |opt| 
			if (opt.eql? option)
				@isfind = 1
				break
			end
		end
		if (@isfind == 1)
			puts "#{name} : #{option}\n\n"
		else
			puts "#{name} : Can't find option <#{option}>"
			print "List Of Option :"
			options.each.with_index { |x, i| print " <#{x}>"}
			puts "\n\n"
			exit
		end
	end
	return option
end

def selectShow(mode, name, option, options = "")
	if (mode.eql? "UI")
		return select(name, options)
	elsif (mode.eql? "CONFIG")
		return show(name, option, options)
	else
		puts "Can't find mode <#{mode}>\n\n"
		exit
	end
end
