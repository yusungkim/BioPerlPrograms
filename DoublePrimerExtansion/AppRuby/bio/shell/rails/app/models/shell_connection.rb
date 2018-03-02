# Model object for rails
class ShellConnection

# include DRb::DRbObservable

  attr_reader :registry

  def puts_remote(str)
    STDOUT.puts(str)
  end

  def initialize
    @connected = false
    @registry = {}
  end

  def [] o_id
    @registry[o_id]
  end

  def []=(name, obj)
    @registry[name] = obj
  end

  def delete(name)
    @registry.delete(name)
  end

end

