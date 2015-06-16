
public class HelloWorld {
	public static void main(String[] args) {

		int len;
		len = args.length;

		System.out.println("Hello, World!!");
		System.out.println(len);

		if (len>0) {
			System.out.print(args[0]);
		    System.out.println(". How are you?");
		
			for (String arg:args) {
				System.out.println(arg);
			}
		}
	}
}
