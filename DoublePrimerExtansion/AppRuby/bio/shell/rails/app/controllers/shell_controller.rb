class DocumentRegistry

  @@registry = []

  def self.[](class_name)
    @@registry.each do |block|
      url = block.call(class_name)
      return url if url
    end
    ""
  end

  def self.register &block
    @@registry.push(block)
  end

end

DocumentRegistry.register() do |class_name|
  if m = /Bio::(.+)/.match(class_name)
    "http://bioruby.org/rdoc/classes/Bio/#{m[1].split('::').join('/')}.html"
  else
    false
  end
end

DocumentRegistry.register do |class_name|
  if m = /Chem::(.+)/.match(class_name)
    "http://chemruby.org/rdoc/classes/Chem/#{m[1].split('::').join('/')}.html"
  else
    false
  end
end

DocumentRegistry.register do |class_name|
  "http://www.ruby-doc.org/core/classes/#{class_name.split('::').join('/')}.html"
end


class ShellController < ApplicationController

  layout 'shell'  #, :except => [:rss_feed, :rss_with_content]

  def index
    setup
  end

  def show
    setup
    @obj = $drb_server[params[:id]]
    if @obj.nil?
      @inheritance = "<h1>Unknown local variable! : #{params[:var]}</h1>"
    else
      @inheritance = get_inheritance @obj.class
      if @obj.respond_to?(:to_html)
        @contents = @obj.to_html
      else
        @contents = "Undefined :to_html"
      end
    end
    @title = params[:id]
  end

=begin
  def history
    setup
    @history = File.read("session/history")
  end
=end

  private

  def setup
    @local_vars = $drb_server.registry
  end
  
  def get_inheritance obj
    #    mods = obj.included_modules - [PP::ObjectMixin, Bio::Shell, WEBrick]
    mods = obj.included_modules
    module_links = mods.collect{|m|
      " [<a href=#{DocumentRegistry[m.to_s]}>#{m.to_s}</a>] "
    }.join("|")

    inherit = []
    loop do
      inherit.push(" [<a href=#{DocumentRegistry[obj.to_s]}>#{obj.to_s}</a>] ")
      break if obj == Object
      obj = obj.superclass
    end
    "<table><tr><th>Inheritance</th><td>" +  inherit.join(" &lt; ") + "</td></tr>" +
      "<tr><th>Mix-in</th><td>" + module_links + "</td</tr></table>"
  end

end
